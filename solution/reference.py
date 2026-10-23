#!/usr/bin/env python3
"""Recover a five-layer nonlinear wall from multi-band thermal records."""

from __future__ import annotations

import argparse
import bisect
import hashlib
import itertools
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


TIME_EPSILON = 1e-7
INTEGRATION_STEP = 0.02
MAX_FIT_ERROR = 0.020
MAX_HEAT_BALANCE_ERROR = 0.5


class ContractError(ValueError):
    """The input or existing ledger violates the frozen contract."""


class DuplicateKeyError(ContractError):
    """A JSON object repeats a key."""


def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKeyError(key)
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    try:
        raw = path.read_bytes().decode("utf-8")
        value = json.loads(raw, object_pairs_hook=pairs_hook)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError("input") from exc
    if not isinstance(value, dict):
        raise ContractError("root")
    return value


def exact_keys(value: dict[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ContractError("keys")


def number(value: Any, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError("number")
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < low or parsed > high:
        raise ContractError("range")
    return parsed


def integer(value: Any, expected: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError("integer")
    if expected is not None and value != expected:
        raise ContractError("integer")
    return value


def label(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise ContractError("label")
    return value


def parse_samples(raw: Any, segment: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, list) or not 48 <= len(raw) <= 120:
        raise ContractError("samples")
    times: list[float] = []
    fluxes: list[float] = []
    observed: list[list[float]] = []
    for item in raw:
        if not isinstance(item, dict):
            raise ContractError("sample")
        exact_keys(item, {"time_s", "heat_flux_w_m2", "sensor_temperatures_k"})
        times.append(number(item["time_s"], 0.0, 30.0))
        fluxes.append(number(item["heat_flux_w_m2"], 1.0, 2000.0))
        sensors = item["sensor_temperatures_k"]
        if not isinstance(sensors, list) or len(sensors) != 5:
            raise ContractError("sensor_count")
        observed.append([number(value, -50.0, 200.0) for value in sensors])
    if abs(times[0] - segment["start_time_s"]) > TIME_EPSILON:
        raise ContractError("segment_origin")
    if abs(times[-1] - segment["end_time_s"]) > TIME_EPSILON:
        raise ContractError("segment_end")
    spacing = times[1] - times[0]
    if spacing <= 0.0:
        raise ContractError("time_order")
    for index in range(1, len(times)):
        if times[index] <= times[index - 1]:
            raise ContractError("time_order")
        if index >= 2 and abs(
            times[index] - times[index - 1] - spacing
        ) > TIME_EPSILON:
            raise ContractError("time_spacing")
    integral = 0.0
    for index in range(1, len(times)):
        integral += 0.5 * (fluxes[index - 1] + fluxes[index]) * (
            times[index] - times[index - 1]
        )
    return {
        "times": times,
        "fluxes": fluxes,
        "observed": observed,
        "spacing": spacing,
        "heat": integral,
    }


def validate(payload: dict[str, Any]) -> dict[str, Any]:
    exact_keys(
        payload,
        {"schema_version", "case_id", "coupons", "segments", "aggregate", "forecast", "limits"},
    )
    integer(payload["schema_version"], 1)
    case_id = label(payload["case_id"])

    limits = payload["limits"]
    if not isinstance(limits, dict):
        raise ContractError("limits")
    exact_keys(
        limits,
        {
            "max_layers",
            "max_coupons",
            "max_segments",
            "max_samples_per_segment",
            "max_forecast_points",
            "max_time_s",
        },
    )
    integer(limits["max_layers"], 5)
    integer(limits["max_coupons"], 5)
    integer(limits["max_segments"], 3)
    integer(limits["max_samples_per_segment"], 120)
    integer(limits["max_forecast_points"], 40)
    number(limits["max_time_s"], 30.0, 30.0)

    raw_coupons = payload["coupons"]
    if not isinstance(raw_coupons, list) or len(raw_coupons) != 5:
        raise ContractError("coupons")
    coupons: list[dict[str, Any]] = []
    coupon_ids: set[str] = set()
    signatures: set[tuple[float, float, float, float]] = set()
    for item in raw_coupons:
        if not isinstance(item, dict):
            raise ContractError("coupon")
        exact_keys(
            item,
            {
                "coupon_id",
                "thickness_m",
                "areal_heat_capacity_j_m2k",
                "k_ref_w_mk",
                "temperature_coefficient_per_k",
            },
        )
        coupon_id = label(item["coupon_id"])
        if coupon_id in coupon_ids:
            raise ContractError("coupon_id")
        coupon_ids.add(coupon_id)
        thickness = number(item["thickness_m"], 0.008, 0.080)
        capacity = number(
            item["areal_heat_capacity_j_m2k"], 400.0, 5000.0
        )
        k_ref = number(item["k_ref_w_mk"], 0.20, 10.0)
        coefficient = number(
            item["temperature_coefficient_per_k"], -0.050, 0.050
        )
        signature = (thickness, capacity, k_ref, coefficient)
        if signature in signatures:
            raise ContractError("duplicate_signature")
        signatures.add(signature)
        coupons.append(
            {
                "coupon_id": coupon_id,
                "thickness_m": thickness,
                "areal_heat_capacity_j_m2k": capacity,
                "k_ref_w_mk": k_ref,
                "temperature_coefficient_per_k": coefficient,
            }
        )

    raw_segments = payload["segments"]
    if not isinstance(raw_segments, list) or len(raw_segments) != 3:
        raise ContractError("segments")
    segments: list[dict[str, Any]] = []
    last_end = -1.0
    total_samples = 0
    for expected_index, item in enumerate(raw_segments, start=1):
        if not isinstance(item, dict):
            raise ContractError("segment")
        exact_keys(
            item,
            {
                "segment_index",
                "frequency_hz",
                "phase_rad",
                "start_time_s",
                "end_time_s",
                "samples",
            },
        )
        integer(item["segment_index"], expected_index)
        frequency = number(item["frequency_hz"], 0.03, 0.50)
        phase = number(item["phase_rad"], 0.0, 2.0 * math.pi)
        start = number(item["start_time_s"], 0.0, 30.0)
        end = number(item["end_time_s"], 0.0, 30.0)
        if end <= start or start <= last_end:
            raise ContractError("segment_time")
        if expected_index == 1 and abs(start) > TIME_EPSILON:
            raise ContractError("segment_origin")
        segment = {
            "segment_index": expected_index,
            "frequency_hz": frequency,
            "phase_rad": phase,
            "start_time_s": start,
            "end_time_s": end,
        }
        parsed_samples = parse_samples(item["samples"], segment)
        segment["data"] = parsed_samples
        segments.append(segment)
        total_samples += len(parsed_samples["times"])
        last_end = end

    aggregate = payload["aggregate"]
    if not isinstance(aggregate, dict):
        raise ContractError("aggregate")
    exact_keys(aggregate, {"total_heat_j_m2", "sample_count"})
    total_heat = number(aggregate["total_heat_j_m2"], 0.0, 1.0e9)
    integer(aggregate["sample_count"], total_samples)
    computed_heat = sum(segment["data"]["heat"] for segment in segments)
    if abs(total_heat - computed_heat) > MAX_HEAT_BALANCE_ERROR:
        raise ContractError("heat_balance")

    forecast = payload["forecast"]
    if not isinstance(forecast, dict):
        raise ContractError("forecast")
    exact_keys(forecast, {"times_s"})
    raw_times = forecast["times_s"]
    if not isinstance(raw_times, list) or not 1 <= len(raw_times) <= 40:
        raise ContractError("forecast_times")
    forecast_times = [number(value, 0.0, 30.0) for value in raw_times]
    if any(
        forecast_times[index] <= forecast_times[index - 1]
        for index in range(1, len(forecast_times))
    ):
        raise ContractError("forecast_order")
    if forecast_times[0] <= last_end:
        raise ContractError("forecast_start")
    if forecast_times[-1] - last_end > 10.0 + TIME_EPSILON:
        raise ContractError("forecast_horizon")

    return {
        "case_id": case_id,
        "coupons": coupons,
        "segments": segments,
        "aggregate_heat": total_heat,
        "forecast_times": forecast_times,
        "final_time": last_end,
    }


def heat_flux_at(segments: list[dict[str, Any]], time_s: float) -> float:
    for segment in segments:
        start = segment["start_time_s"]
        end = segment["end_time_s"]
        if time_s < start or time_s > end:
            continue
        data = segment["data"]
        times = data["times"]
        index = bisect.bisect_right(times, time_s) - 1
        index = max(0, min(index, len(times) - 1))
        if index == len(times) - 1:
            return data["fluxes"][index]
        left = times[index]
        right = times[index + 1]
        ratio = 0.0 if right == left else (time_s - left) / (right - left)
        return data["fluxes"][index] + ratio * (
            data["fluxes"][index + 1] - data["fluxes"][index]
        )
    return 0.0


def conductance(
    coupon: dict[str, Any], left_temperature: float, right_temperature: float
) -> float:
    average = 0.5 * (left_temperature + right_temperature)
    value = coupon["k_ref_w_mk"] * (
        1.0 + coupon["temperature_coefficient_per_k"] * average
    )
    if value <= 0.0 or not math.isfinite(value):
        raise ContractError("conductivity")
    return value / coupon["thickness_m"]


def derivative(
    state: list[float],
    ordered_coupons: list[dict[str, Any]],
    flux: float,
) -> list[float]:
    conductances: list[float] = []
    for index, coupon in enumerate(ordered_coupons):
        left = 0.0 if index == 0 else state[index - 1]
        conductances.append(conductance(coupon, left, state[index]))
    result: list[float] = []
    for index, coupon in enumerate(ordered_coupons):
        capacity = coupon["areal_heat_capacity_j_m2k"]
        left_temperature = 0.0 if index == 0 else state[index - 1]
        left_flow = conductances[index] * (left_temperature - state[index])
        if index == len(ordered_coupons) - 1:
            right_flow = flux
        else:
            right_flow = conductances[index + 1] * (
                state[index + 1] - state[index]
            )
        result.append((left_flow + right_flow) / capacity)
    return result


def rk4_step(
    state: list[float],
    ordered_coupons: list[dict[str, Any]],
    time_s: float,
    step: float,
    flux_function: Callable[[float], float],
) -> list[float]:
    k1 = derivative(state, ordered_coupons, flux_function(time_s))
    state2 = [
        state[index] + 0.5 * step * k1[index] for index in range(len(state))
    ]
    k2 = derivative(
        state2, ordered_coupons, flux_function(time_s + 0.5 * step)
    )
    state3 = [
        state[index] + 0.5 * step * k2[index] for index in range(len(state))
    ]
    k3 = derivative(
        state3, ordered_coupons, flux_function(time_s + 0.5 * step)
    )
    state4 = [
        state[index] + step * k3[index] for index in range(len(state))
    ]
    k4 = derivative(
        state4, ordered_coupons, flux_function(time_s + step)
    )
    return [
        state[index]
        + step
        * (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index])
        / 6.0
        for index in range(len(state))
    ]


def advance(
    state: list[float],
    ordered_coupons: list[dict[str, Any]],
    time_s: float,
    target_time: float,
    flux_function: Callable[[float], float],
) -> tuple[list[float], float]:
    current_state = list(state)
    current_time = time_s
    while current_time < target_time - TIME_EPSILON:
        step = min(INTEGRATION_STEP, target_time - current_time)
        current_state = rk4_step(
            current_state, ordered_coupons, current_time, step, flux_function
        )
        current_time += step
    return current_state, target_time


def simulate_observations(
    parsed: dict[str, Any], order: tuple[int, ...]
) -> tuple[list[list[float]], list[float], float]:
    coupons = parsed["coupons"]
    ordered_coupons = [coupons[index] for index in order]
    flux_function = lambda time_s: heat_flux_at(parsed["segments"], time_s)
    state = [0.0] * 5
    current_time = 0.0
    observations: list[list[float]] = []
    for segment in parsed["segments"]:
        state, current_time = advance(
            state,
            ordered_coupons,
            current_time,
            segment["start_time_s"],
            flux_function,
        )
        for sample_time in segment["data"]["times"]:
            state, current_time = advance(
                state,
                ordered_coupons,
                current_time,
                sample_time,
                flux_function,
            )
            observations.append(list(state))
        state, current_time = advance(
            state,
            ordered_coupons,
            current_time,
            segment["end_time_s"],
            flux_function,
        )
    return observations, state, current_time


def observed_matrix(parsed: dict[str, Any]) -> list[list[float]]:
    rows: list[list[float]] = []
    for segment in parsed["segments"]:
        rows.extend(segment["data"]["observed"])
    return rows


def fit_calibration(
    predicted: list[list[float]], observed: list[list[float]]
) -> tuple[list[dict[str, float]], float]:
    calibrations: list[dict[str, float]] = []
    maximum_error = 0.0
    for sensor in range(5):
        temperatures = [row[sensor] for row in predicted]
        readings = [row[sensor] for row in observed]
        count = len(temperatures)
        mean_x = sum(temperatures) / count
        mean_y = sum(readings) / count
        denominator = sum(
            (temperature - mean_x) ** 2 for temperature in temperatures
        )
        if denominator <= 1e-14:
            raise ContractError("sensor_variance")
        gain = (
            sum(
                (temperature - mean_x) * (reading - mean_y)
                for temperature, reading in zip(temperatures, readings)
            )
            / denominator
        )
        offset = mean_y - gain * mean_x
        if not 0.70 <= gain <= 1.30 or not -0.30 <= offset <= 0.30:
            raise ContractError("calibration_range")
        for temperature, reading in zip(temperatures, readings):
            maximum_error = max(
                maximum_error, abs(gain * temperature + offset - reading)
            )
        calibrations.append({"gain": gain, "offset_k": offset})
    return calibrations, maximum_error


def recover(parsed: dict[str, Any]) -> tuple[
    tuple[int, ...],
    list[dict[str, float]],
    list[list[float]],
    list[float],
    float,
    float,
    float,
]:
    observed = observed_matrix(parsed)
    candidates: list[
        tuple[float, float, tuple[int, ...], list[dict[str, float]], list[float]]
    ] = []
    for order in itertools.permutations(range(5)):
        predicted, final_state, final_time = simulate_observations(parsed, order)
        try:
            calibration, maximum_error = fit_calibration(predicted, observed)
        except ContractError:
            continue
        sse = 0.0
        for predicted_row, observed_row in zip(predicted, observed):
            for sensor in range(5):
                residual = (
                    calibration[sensor]["gain"] * predicted_row[sensor]
                    + calibration[sensor]["offset_k"]
                    - observed_row[sensor]
                )
                sse += residual * residual
        candidates.append(
            (sse, maximum_error, order, calibration, final_state)
        )
    if not candidates:
        raise ContractError("no_candidate")
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    best = candidates[0]
    if best[1] > MAX_FIT_ERROR:
        raise ContractError("fit")
    if len(candidates) > 1:
        second = candidates[1]
        required_margin = max(1.0e-8, 8.0 * best[0])
        if second[0] - best[0] <= required_margin:
            raise ContractError("ambiguous")

    order = best[2]
    calibration = best[3]
    final_state = best[4]
    flux_function = lambda time_s: 0.0
    forecast_rows: list[list[float]] = []
    state = list(final_state)
    current_time = parsed["final_time"]
    for forecast_time in parsed["forecast_times"]:
        state, current_time = advance(
            state,
            [parsed["coupons"][index] for index in order],
            current_time,
            forecast_time,
            flux_function,
        )
        forecast_rows.append(
            [
                calibration[sensor]["gain"] * state[sensor]
                + calibration[sensor]["offset_k"]
                for sensor in range(5)
            ]
        )
    total_heat = sum(
        segment["data"]["heat"] for segment in parsed["segments"]
    )
    heat_error = parsed["aggregate_heat"] - total_heat
    return (
        order,
        calibration,
        forecast_rows,
        forecast_times(parsed),
        total_heat,
        heat_error,
        best[1],
    )


def forecast_times(parsed: dict[str, Any]) -> list[float]:
    return parsed["forecast_times"]


def build_result(parsed: dict[str, Any]) -> dict[str, Any]:
    (
        order,
        calibration,
        forecast_rows,
        times,
        total_heat,
        heat_error,
        maximum_error,
    ) = recover(parsed)
    coupons = parsed["coupons"]
    layer_order = [
        {
            "slot_index": slot,
            "coupon_id": coupons[index]["coupon_id"],
        }
        for slot, index in enumerate(order, start=1)
    ]
    conductivity_curves = [
        {
            "slot_index": slot,
            "k_ref_w_mk": round(coupons[index]["k_ref_w_mk"], 9),
            "temperature_coefficient_per_k": round(
                coupons[index]["temperature_coefficient_per_k"], 9
            ),
        }
        for slot, index in enumerate(order, start=1)
    ]
    sensor_calibration = [
        {
            "sensor_index": sensor,
            "gain": round(item["gain"], 9),
            "offset_k": round(item["offset_k"], 9),
        }
        for sensor, item in enumerate(calibration, start=1)
    ]
    forecast = [
        {
            "sensor_index": sensor,
            "samples": [
                {
                    "time_s": round(times[index], 9),
                    "temperature_rise_k": round(row[sensor - 1], 9),
                }
                for index, row in enumerate(forecast_rows)
            ],
        }
        for sensor in range(1, 6)
    ]
    return {
        "schema_version": 1,
        "case_id": parsed["case_id"],
        "layer_order": layer_order,
        "conductivity_curves": conductivity_curves,
        "sensor_calibration": sensor_calibration,
        "forecast": forecast,
        "fit": {
            "max_abs_temperature_error_k": round(maximum_error, 9),
            "total_heat_j_m2": round(total_heat, 6),
            "heat_balance_error_j_m2": round(heat_error, 6),
        },
    }


def canonical_bytes(value: dict[str, Any]) -> bytes:
    text = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return (text + "\n").encode("utf-8")


def validate_ledger(path: Path) -> bytes:
    if not path.exists():
        return b""
    if not path.is_file():
        raise ContractError("ledger_type")
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        raise ContractError("ledger_newline")
    for line in raw.splitlines():
        try:
            record = json.loads(line.decode("utf-8"), object_pairs_hook=pairs_hook)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise ContractError("ledger_json") from exc
        if not isinstance(record, dict):
            raise ContractError("ledger_record")
        exact_keys(
            record,
            {
                "case_id",
                "event",
                "input_sha256",
                "result_sha256",
                "schema_version",
            },
        )
        if record["event"] != "completed" or record["schema_version"] != 1:
            raise ContractError("ledger_record")
    return raw


def stage_bytes(path: Path, data: bytes, mode: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary)
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(temporary_path, mode)
    return temporary_path


def install_outputs(
    output_path: Path,
    output_bytes: bytes,
    ledger_path: Path,
    ledger_bytes: bytes | None,
) -> None:
    if output_path == ledger_path:
        raise ContractError("paths")
    if output_path.exists() and not output_path.is_file():
        raise ContractError("output_type")
    if ledger_path.exists() and not ledger_path.is_file():
        raise ContractError("ledger_type")
    old_output = output_path.read_bytes() if output_path.exists() else None
    output_temp = stage_bytes(output_path, output_bytes, 0o644)
    ledger_temp: Path | None = None
    try:
        if ledger_bytes is not None:
            ledger_temp = stage_bytes(ledger_path, ledger_bytes, 0o600)
        os.replace(output_temp, output_path)
        if ledger_temp is not None:
            try:
                os.replace(ledger_temp, ledger_path)
            except Exception:
                if old_output is None:
                    output_path.unlink(missing_ok=True)
                else:
                    restore = stage_bytes(output_path, old_output, 0o644)
                    os.replace(restore, output_path)
                raise
    finally:
        output_temp.unlink(missing_ok=True)
        if ledger_temp is not None:
            ledger_temp.unlink(missing_ok=True)


def execute(args: argparse.Namespace) -> None:
    input_path = Path(args.input)
    output_path = Path(args.output)
    ledger_path = Path(args.ledger)
    parsed = validate(load_json(input_path))
    result = build_result(parsed)
    input_bytes = input_path.read_bytes()
    result_bytes = canonical_bytes(result)
    record = {
        "schema_version": 1,
        "event": "completed",
        "case_id": parsed["case_id"],
        "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "result_sha256": hashlib.sha256(result_bytes).hexdigest(),
    }
    line = canonical_bytes(record)
    existing = validate_ledger(ledger_path)
    if line in existing.splitlines(keepends=True):
        planned = existing
        ledger_change: bytes | None = None
    else:
        planned = existing + line
        ledger_change = planned
    install_outputs(output_path, result_bytes, ledger_path, ledger_change)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ledger", required=True)
    return parser.parse_args(argv)


def main() -> int:
    try:
        execute(parse_args(sys.argv[1:]))
    except Exception:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
