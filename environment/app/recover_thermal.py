#!/usr/bin/env python3
"""Previous wall reporter retained for the furnace records."""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


STEP = 0.02
EPSILON = 1e-7


class InputProblem(ValueError):
    """The thermal record is not usable."""


class DuplicateKey(InputProblem):
    """A JSON object repeats a key."""


def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKey(key)
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_bytes().decode("utf-8"), object_pairs_hook=pairs_hook
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise InputProblem("read") from exc
    if not isinstance(value, dict):
        raise InputProblem("object")
    return value


def exact_keys(value: dict[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise InputProblem("keys")


def numeric(value: Any, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputProblem("number")
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < low or parsed > high:
        raise InputProblem("range")
    return parsed


def whole(value: Any, expected: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InputProblem("integer")
    if expected is not None and value != expected:
        raise InputProblem("integer")
    return value


def label(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise InputProblem("label")
    return value


def parse_samples(raw: Any, segment: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, list) or not 48 <= len(raw) <= 120:
        raise InputProblem("samples")
    times: list[float] = []
    fluxes: list[float] = []
    observed: list[list[float]] = []
    for item in raw:
        if not isinstance(item, dict):
            raise InputProblem("sample")
        exact_keys(item, {"time_s", "heat_flux_w_m2", "sensor_temperatures_k"})
        times.append(numeric(item["time_s"], 0.0, 30.0))
        fluxes.append(numeric(item["heat_flux_w_m2"], 1.0, 2000.0))
        sensors = item["sensor_temperatures_k"]
        if not isinstance(sensors, list) or len(sensors) != 5:
            raise InputProblem("sensors")
        observed.append([numeric(value, -50.0, 200.0) for value in sensors])
    if abs(times[0] - segment["start_time_s"]) > EPSILON:
        raise InputProblem("start")
    if abs(times[-1] - segment["end_time_s"]) > EPSILON:
        raise InputProblem("end")
    spacing = times[1] - times[0]
    if any(times[index] <= times[index - 1] for index in range(1, len(times))):
        raise InputProblem("time_order")
    if any(
        abs(times[index] - times[index - 1] - spacing) > EPSILON
        for index in range(2, len(times))
    ):
        raise InputProblem("spacing")
    heat = sum(
        0.5
        * (fluxes[index - 1] + fluxes[index])
        * (times[index] - times[index - 1])
        for index in range(1, len(times))
    )
    return {
        "times": times,
        "fluxes": fluxes,
        "observed": observed,
        "heat": heat,
    }


def validate(payload: dict[str, Any]) -> dict[str, Any]:
    exact_keys(
        payload,
        {"schema_version", "case_id", "coupons", "segments", "aggregate", "forecast", "limits"},
    )
    whole(payload["schema_version"], 1)
    case_id = label(payload["case_id"])
    limits = payload["limits"]
    if not isinstance(limits, dict):
        raise InputProblem("limits")
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
    whole(limits["max_layers"], 5)
    whole(limits["max_coupons"], 5)
    whole(limits["max_segments"], 3)
    whole(limits["max_samples_per_segment"], 120)
    whole(limits["max_forecast_points"], 40)
    numeric(limits["max_time_s"], 30.0, 30.0)

    raw_coupons = payload["coupons"]
    if not isinstance(raw_coupons, list) or len(raw_coupons) != 5:
        raise InputProblem("coupons")
    coupons: list[dict[str, Any]] = []
    coupon_ids: set[str] = set()
    signatures: set[tuple[float, float, float, float]] = set()
    for item in raw_coupons:
        if not isinstance(item, dict):
            raise InputProblem("coupon")
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
            raise InputProblem("coupon_id")
        coupon_ids.add(coupon_id)
        thickness = numeric(item["thickness_m"], 0.008, 0.080)
        capacity = numeric(
            item["areal_heat_capacity_j_m2k"], 400.0, 5000.0
        )
        k_ref = numeric(item["k_ref_w_mk"], 0.20, 10.0)
        coefficient = numeric(
            item["temperature_coefficient_per_k"], -0.050, 0.050
        )
        signature = (thickness, capacity, k_ref, coefficient)
        if signature in signatures:
            raise InputProblem("signature")
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
        raise InputProblem("segments")
    segments: list[dict[str, Any]] = []
    last_end = -1.0
    total_samples = 0
    for expected_index, item in enumerate(raw_segments, start=1):
        if not isinstance(item, dict):
            raise InputProblem("segment")
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
        whole(item["segment_index"], expected_index)
        frequency = numeric(item["frequency_hz"], 0.03, 0.50)
        phase = numeric(item["phase_rad"], 0.0, 2.0 * math.pi)
        start = numeric(item["start_time_s"], 0.0, 30.0)
        end = numeric(item["end_time_s"], 0.0, 30.0)
        if end <= start or start <= last_end:
            raise InputProblem("segment_time")
        if expected_index == 1 and abs(start) > EPSILON:
            raise InputProblem("origin")
        segment = {
            "segment_index": expected_index,
            "frequency_hz": frequency,
            "phase_rad": phase,
            "start_time_s": start,
            "end_time_s": end,
        }
        segment["data"] = parse_samples(item["samples"], segment)
        segments.append(segment)
        total_samples += len(segment["data"]["times"])
        last_end = end

    aggregate = payload["aggregate"]
    if not isinstance(aggregate, dict):
        raise InputProblem("aggregate")
    exact_keys(aggregate, {"total_heat_j_m2", "sample_count"})
    total_heat = numeric(aggregate["total_heat_j_m2"], 0.0, 1.0e9)
    whole(aggregate["sample_count"], total_samples)
    computed_heat = sum(segment["data"]["heat"] for segment in segments)
    if abs(total_heat - computed_heat) > 0.5:
        raise InputProblem("heat")

    forecast = payload["forecast"]
    if not isinstance(forecast, dict):
        raise InputProblem("forecast")
    exact_keys(forecast, {"times_s"})
    raw_times = forecast["times_s"]
    if not isinstance(raw_times, list) or not 1 <= len(raw_times) <= 40:
        raise InputProblem("forecast_times")
    forecast_times = [numeric(value, 0.0, 30.0) for value in raw_times]
    if any(
        forecast_times[index] <= forecast_times[index - 1]
        for index in range(1, len(forecast_times))
    ):
        raise InputProblem("forecast_order")
    if forecast_times[0] <= last_end:
        raise InputProblem("forecast_start")
    if forecast_times[-1] - last_end > 10.0 + EPSILON:
        raise InputProblem("forecast_horizon")
    return {
        "case_id": case_id,
        "coupons": coupons,
        "segments": segments,
        "forecast_times": forecast_times,
        "final_time": last_end,
        "aggregate_heat": total_heat,
    }


def flux_at(segments: list[dict[str, Any]], time_s: float) -> float:
    for segment in segments:
        if time_s < segment["start_time_s"] or time_s > segment["end_time_s"]:
            continue
        times = segment["data"]["times"]
        values = segment["data"]["fluxes"]
        index = max(0, min(bisect.bisect_right(times, time_s) - 1, len(times) - 1))
        if index == len(times) - 1:
            return values[index]
        ratio = (time_s - times[index]) / (times[index + 1] - times[index])
        return values[index] + ratio * (values[index + 1] - values[index])
    return 0.0


def derivative(
    state: list[float], coupons: list[dict[str, Any]], flux: float
) -> list[float]:
    conductance = [
        coupon["k_ref_w_mk"] / coupon["thickness_m"] for coupon in coupons
    ]
    result: list[float] = []
    for index, coupon in enumerate(coupons):
        left_temperature = 0.0 if index == 0 else state[index - 1]
        left_flow = conductance[index] * (left_temperature - state[index])
        right_flow = (
            flux
            if index == 4
            else conductance[index + 1] * (state[index + 1] - state[index])
        )
        result.append(
            (left_flow + right_flow) / coupon["areal_heat_capacity_j_m2k"]
        )
    return result


def advance(
    state: list[float],
    coupons: list[dict[str, Any]],
    time_s: float,
    target: float,
    flux_function: Callable[[float], float],
) -> tuple[list[float], float]:
    current = list(state)
    while time_s < target - EPSILON:
        step = min(STEP, target - time_s)
        k1 = derivative(current, coupons, flux_function(time_s))
        second = [current[i] + 0.5 * step * k1[i] for i in range(5)]
        k2 = derivative(second, coupons, flux_function(time_s + 0.5 * step))
        third = [current[i] + 0.5 * step * k2[i] for i in range(5)]
        k3 = derivative(third, coupons, flux_function(time_s + 0.5 * step))
        fourth = [current[i] + step * k3[i] for i in range(5)]
        k4 = derivative(fourth, coupons, flux_function(time_s + step))
        current = [
            current[i]
            + step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            for i in range(5)
        ]
        time_s += step
    return current, target


def simulate(
    parsed: dict[str, Any], coupons: list[dict[str, Any]]
) -> tuple[list[list[float]], list[float]]:
    flux_function = lambda time_s: flux_at(parsed["segments"], time_s)
    state = [0.0] * 5
    current_time = 0.0
    rows: list[list[float]] = []
    for segment in parsed["segments"]:
        state, current_time = advance(
            state,
            coupons,
            current_time,
            segment["start_time_s"],
            flux_function,
        )
        for sample_time in segment["data"]["times"]:
            state, current_time = advance(
                state, coupons, current_time, sample_time, flux_function
            )
            rows.append(list(state))
        state, current_time = advance(
            state,
            coupons,
            current_time,
            segment["end_time_s"],
            flux_function,
        )
    return rows, state


def fit_calibration(
    predicted: list[list[float]], observed: list[list[float]]
) -> list[dict[str, float]]:
    result: list[dict[str, float]] = []
    for sensor in range(5):
        x = [row[sensor] for row in predicted]
        y = [row[sensor] for row in observed]
        mean_x = sum(x) / len(x)
        mean_y = sum(y) / len(y)
        denominator = sum((value - mean_x) ** 2 for value in x)
        if denominator <= 1e-12:
            result.append({"gain": 1.0, "offset_k": 0.0})
            continue
        gain = sum(
            (left - mean_x) * (right - mean_y) for left, right in zip(x, y)
        ) / denominator
        offset = mean_y - gain * mean_x
        result.append({"gain": gain, "offset_k": offset})
    return result


def solve(parsed: dict[str, Any]) -> dict[str, Any]:
    order = tuple(
        sorted(
            range(5),
            key=lambda index: parsed["coupons"][index]["thickness_m"]
            / parsed["coupons"][index]["k_ref_w_mk"],
        )
    )
    coupons = [parsed["coupons"][index] for index in order]
    predicted, state = simulate(parsed, coupons)
    observed = [
        row
        for segment in parsed["segments"]
        for row in segment["data"]["observed"]
    ]
    calibration = fit_calibration(predicted, observed)
    error = max(
        abs(
            calibration[sensor]["gain"] * predicted[index][sensor]
            + calibration[sensor]["offset_k"]
            - observed[index][sensor]
        )
        for index in range(len(observed))
        for sensor in range(5)
    )
    forecast_rows: list[list[float]] = []
    current_time = parsed["final_time"]
    for forecast_time in parsed["forecast_times"]:
        state, current_time = advance(
            state, coupons, current_time, forecast_time, lambda _: 0.0
        )
        forecast_rows.append(
            [
                calibration[sensor]["gain"] * state[sensor]
                + calibration[sensor]["offset_k"]
                for sensor in range(5)
            ]
        )
    return {
        "schema_version": 1,
        "case_id": parsed["case_id"],
        "layer_order": [
            {
                "slot_index": slot,
                "coupon_id": parsed["coupons"][index]["coupon_id"],
            }
            for slot, index in enumerate(order, start=1)
        ],
        "conductivity_curves": [
            {
                "slot_index": slot,
                "k_ref_w_mk": round(coupons[slot - 1]["k_ref_w_mk"], 9),
                "temperature_coefficient_per_k": round(
                    coupons[slot - 1]["temperature_coefficient_per_k"], 9
                ),
            }
            for slot in range(1, 6)
        ],
        "sensor_calibration": [
            {
                "sensor_index": sensor,
                "gain": round(item["gain"], 9),
                "offset_k": round(item["offset_k"], 9),
            }
            for sensor, item in enumerate(calibration, start=1)
        ],
        "forecast": [
            {
                "sensor_index": sensor,
                "samples": [
                    {
                        "time_s": round(parsed["forecast_times"][index], 9),
                        "temperature_rise_k": round(row[sensor - 1], 9),
                    }
                    for index, row in enumerate(forecast_rows)
                ],
            }
            for sensor in range(1, 6)
        ],
        "fit": {
            "max_abs_temperature_error_k": round(error, 9),
            "total_heat_j_m2": round(
                sum(
                    segment["data"]["heat"]
                    for segment in parsed["segments"]
                ),
                6,
            ),
            "heat_balance_error_j_m2": round(
                parsed["aggregate_heat"]
                - sum(
                    segment["data"]["heat"]
                    for segment in parsed["segments"]
                ),
                6,
            ),
        },
    }


def compact(value: dict[str, Any]) -> bytes:
    text = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return (text + "\n").encode("utf-8")


def validate_ledger(path: Path) -> bytes:
    if not path.exists():
        return b""
    if not path.is_file():
        raise InputProblem("ledger_type")
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        raise InputProblem("ledger_newline")
    for line in raw.splitlines():
        try:
            record = json.loads(line.decode("utf-8"), object_pairs_hook=pairs_hook)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise InputProblem("ledger_json") from exc
        if not isinstance(record, dict):
            raise InputProblem("ledger")
        exact_keys(
            record,
            {"case_id", "event", "input_sha256", "result_sha256", "schema_version"},
        )
    return raw


def stage(path: Path, data: bytes, mode: int) -> Path:
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


def install(
    output_path: Path,
    output_bytes: bytes,
    ledger_path: Path,
    ledger_bytes: bytes | None,
) -> None:
    if output_path.exists() and not output_path.is_file():
        raise InputProblem("output_type")
    if ledger_path.exists() and not ledger_path.is_file():
        raise InputProblem("ledger_type")
    old_output = output_path.read_bytes() if output_path.exists() else None
    output_temp = stage(output_path, output_bytes, 0o644)
    ledger_temp: Path | None = None
    try:
        if ledger_bytes is not None:
            ledger_temp = stage(ledger_path, ledger_bytes, 0o600)
        os.replace(output_temp, output_path)
        if ledger_temp is not None:
            os.replace(ledger_temp, ledger_path)
    except Exception:
        if old_output is None:
            output_path.unlink(missing_ok=True)
        else:
            restore = stage(output_path, old_output, 0o644)
            os.replace(restore, output_path)
        raise
    finally:
        output_temp.unlink(missing_ok=True)
        if ledger_temp is not None:
            ledger_temp.unlink(missing_ok=True)


def execute(args: argparse.Namespace) -> None:
    input_path = Path(args.input)
    parsed = validate(load_json(input_path))
    result = solve(parsed)
    result_bytes = compact(result)
    record = {
        "schema_version": 1,
        "event": "completed",
        "case_id": parsed["case_id"],
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "result_sha256": hashlib.sha256(result_bytes).hexdigest(),
    }
    line = compact(record)
    ledger_path = Path(args.ledger)
    existing = validate_ledger(ledger_path)
    ledger_change = None if line in existing.splitlines(keepends=True) else existing + line
    install(Path(args.output), result_bytes, ledger_path, ledger_change)


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
