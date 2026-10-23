from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path

import pytest


ENTRY = os.environ.get("YEL006_ENTRY", "/app/recover_thermal.py")
PARAMETERS = [
    (0.024, 900.0, 1.8, 0.018),
    (0.032, 1200.0, 2.4, -0.012),
    (0.018, 650.0, 1.2, 0.030),
    (0.041, 1500.0, 3.0, -0.018),
    (0.027, 1050.0, 1.5, 0.026),
]
TRUE_ORDER = (4, 1, 3, 0, 2)
TRUE_GAINS = [0.87, 1.13, 0.93, 1.08, 0.96]
TRUE_OFFSETS = [-0.14, 0.11, 0.06, -0.09, 0.18]
SEGMENT_DEFINITIONS = [
    (0.0, 4.2, 70, 0.11, 0.2, 420.0),
    (4.5, 9.4, 90, 0.23, 1.1, 560.0),
    (9.8, 15.2, 110, 0.37, 2.0, 480.0),
]
PRIMARY_FORECAST_TIMES = [16.0, 17.5, 19.0, 21.0, 23.5]
STEP = 0.02


def _coupons() -> list[dict]:
    return [
        {
            "coupon_id": f"C{index:03d}",
            "thickness_m": thickness,
            "areal_heat_capacity_j_m2k": capacity,
            "k_ref_w_mk": k_ref,
            "temperature_coefficient_per_k": coefficient,
        }
        for index, (thickness, capacity, k_ref, coefficient) in enumerate(
            PARAMETERS, start=1
        )
    ]


def _fluxes(
    start: float,
    count: int,
    end: float,
    frequency: float,
    phase: float,
    amplitude: float,
) -> tuple[list[float], list[float]]:
    times = [start + (end - start) * index / (count - 1) for index in range(count)]
    values = [
        amplitude
        * (
            0.72
            + 0.20 * math.sin(2.0 * math.pi * frequency * (time_s - start) + phase)
            + 0.08
            * math.sin(4.0 * math.pi * frequency * (time_s - start) + 0.5 * phase)
        )
        for time_s in times
    ]
    return times, values


def _segments(
    sample_counts: tuple[int, int, int] = (70, 90, 110),
) -> tuple[list[dict], float]:
    segments: list[dict] = []
    total_heat = 0.0
    for index, definition in enumerate(SEGMENT_DEFINITIONS):
        start, end, _, frequency, phase, amplitude = definition
        times, fluxes = _fluxes(
            start, sample_counts[index], end, frequency, phase, amplitude
        )
        total_heat += sum(
            0.5
            * (fluxes[position - 1] + fluxes[position])
            * (times[position] - times[position - 1])
            for position in range(1, len(times))
        )
        segments.append(
            {
                "segment_index": index + 1,
                "frequency_hz": frequency,
                "phase_rad": phase,
                "start_time_s": start,
                "end_time_s": end,
                "samples": [
                    {
                        "time_s": round(time_s, 9),
                        "heat_flux_w_m2": round(flux, 9),
                        "sensor_temperatures_k": [0.0] * 5,
                    }
                    for time_s, flux in zip(times, fluxes)
                ],
            }
        )
    return segments, total_heat


def _flux_at(segments: list[dict], time_s: float) -> float:
    for segment in segments:
        if time_s < segment["start_time_s"] or time_s > segment["end_time_s"]:
            continue
        samples = segment["samples"]
        times = [sample["time_s"] for sample in samples]
        for index in range(len(times) - 1):
            if times[index] <= time_s <= times[index + 1]:
                ratio = (time_s - times[index]) / (times[index + 1] - times[index])
                return samples[index]["heat_flux_w_m2"] + ratio * (
                    samples[index + 1]["heat_flux_w_m2"]
                    - samples[index]["heat_flux_w_m2"]
                )
        return samples[-1]["heat_flux_w_m2"]
    return 0.0


def _derivative(
    state: list[float],
    coupons: list[dict],
    flux: float,
) -> list[float]:
    conductance = []
    for index, coupon in enumerate(coupons):
        left = 0.0 if index == 0 else state[index - 1]
        average = 0.5 * (left + state[index])
        coefficient = coupon["temperature_coefficient_per_k"]
        k_value = coupon["k_ref_w_mk"] * (1.0 + coefficient * average)
        conductance.append(k_value / coupon["thickness_m"])
    result = []
    for index, coupon in enumerate(coupons):
        left = 0.0 if index == 0 else state[index - 1]
        left_flow = conductance[index] * (left - state[index])
        right_flow = (
            flux
            if index == 4
            else conductance[index + 1] * (state[index + 1] - state[index])
        )
        result.append(
            (left_flow + right_flow)
            / coupon["areal_heat_capacity_j_m2k"]
        )
    return result


def _advance(
    state: list[float],
    coupons: list[dict],
    time_s: float,
    target: float,
    segments: list[dict],
) -> tuple[list[float], float]:
    current = list(state)
    while time_s < target - 1e-9:
        step = min(STEP, target - time_s)
        flux1 = _flux_at(segments, time_s)
        k1 = _derivative(current, coupons, flux1)
        second = [current[index] + 0.5 * step * k1[index] for index in range(5)]
        k2 = _derivative(second, coupons, _flux_at(segments, time_s + 0.5 * step))
        third = [current[index] + 0.5 * step * k2[index] for index in range(5)]
        k3 = _derivative(third, coupons, _flux_at(segments, time_s + 0.5 * step))
        fourth = [current[index] + step * k3[index] for index in range(5)]
        k4 = _derivative(fourth, coupons, _flux_at(segments, time_s + step))
        current = [
            current[index]
            + step
            * (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index])
            / 6.0
            for index in range(5)
        ]
        time_s += step
    return current, target


def _simulate(
    coupons: list[dict],
    segments: list[dict],
    order: tuple[int, ...],
) -> tuple[list[list[float]], list[float], float]:
    ordered = [coupons[index] for index in order]
    state = [0.0] * 5
    current_time = 0.0
    rows: list[list[float]] = []
    for segment in segments:
        state, current_time = _advance(
            state,
            ordered,
            current_time,
            segment["start_time_s"],
            segments,
        )
        for sample in segment["samples"]:
            state, current_time = _advance(
                state,
                ordered,
                current_time,
                sample["time_s"],
                segments,
            )
            rows.append(list(state))
        state, current_time = _advance(
            state,
            ordered,
            current_time,
            segment["end_time_s"],
            segments,
        )
    return rows, state, current_time


def _payload(
    case_id: str = "wall-nonlinear-0611",
    sample_counts: tuple[int, int, int] = (70, 90, 110),
    forecast_times: list[float] | None = None,
) -> dict:
    coupons = _coupons()
    segments, total_heat = _segments(sample_counts)
    rows, final_state, final_time = _simulate(coupons, segments, TRUE_ORDER)
    position = 0
    for segment in segments:
        for sample in segment["samples"]:
            sample["sensor_temperatures_k"] = [
                round(TRUE_GAINS[sensor] * rows[position][sensor] + TRUE_OFFSETS[sensor], 9)
                for sensor in range(5)
            ]
            position += 1
    times = PRIMARY_FORECAST_TIMES if forecast_times is None else forecast_times
    return {
        "schema_version": 1,
        "case_id": case_id,
        "coupons": coupons,
        "segments": segments,
        "aggregate": {
            "total_heat_j_m2": round(total_heat, 6),
            "sample_count": sum(len(segment["samples"]) for segment in segments),
        },
        "forecast": {"times_s": times},
        "limits": {
            "max_layers": 5,
            "max_coupons": 5,
            "max_segments": 3,
            "max_samples_per_segment": 120,
            "max_forecast_points": 40,
            "max_time_s": 30.0,
        },
    }


def _expected_forecast(payload: dict) -> list[list[float]]:
    observations, state, current_time = _simulate(
        payload["coupons"], payload["segments"], TRUE_ORDER
    )
    del observations
    rows = []
    ordered = [payload["coupons"][index] for index in TRUE_ORDER]
    for forecast_time in payload["forecast"]["times_s"]:
        state, current_time = _advance(
            state, ordered, current_time, forecast_time, []
        )
        rows.append(
            [
                TRUE_GAINS[sensor] * state[sensor] + TRUE_OFFSETS[sensor]
                for sensor in range(5)
            ]
        )
    return rows


def _invoke(
    directory: Path,
    payload: dict,
    existing_output: bytes | None = None,
    existing_ledger: bytes | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path, Path, Path]:
    input_path = directory / "wall-input.json"
    output_path = directory / "wall-result.json"
    ledger_path = directory / "wall-ledger.jsonl"
    input_path.write_text(
        json.dumps(payload, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    if existing_output is not None:
        output_path.write_bytes(existing_output)
    if existing_ledger is not None:
        ledger_path.write_bytes(existing_ledger)
    completed = subprocess.run(
        [
            "python3",
            ENTRY,
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ledger",
            str(ledger_path),
        ],
        capture_output=True,
        text=True,
        timeout=45,
        env=os.environ.copy(),
    )
    return completed, output_path, ledger_path, input_path


def _raw_invoke(
    directory: Path, raw: bytes
) -> tuple[subprocess.CompletedProcess[str], Path, Path]:
    input_path = directory / "raw-input.json"
    output_path = directory / "raw-result.json"
    ledger_path = directory / "raw-ledger.jsonl"
    input_path.write_bytes(raw)
    completed = subprocess.run(
        [
            "python3",
            ENTRY,
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ledger",
            str(ledger_path),
        ],
        capture_output=True,
        text=True,
        timeout=45,
        env=os.environ.copy(),
    )
    return completed, output_path, ledger_path


@pytest.fixture(scope="module")
def solved(tmp_path_factory):
    directory = tmp_path_factory.mktemp("primary-wall")
    payload = _payload()
    completed, output_path, ledger_path, input_path = _invoke(directory, payload)
    raw = output_path.read_bytes() if output_path.exists() else b""
    result = json.loads(raw) if raw else {}
    return {
        "payload": payload,
        "completed": completed,
        "result": result,
        "raw": raw,
        "ledger": ledger_path,
        "input_path": input_path,
        "directory": directory,
    }


def test_recovers_hidden_layer_order(solved):
    """题面：「装配序列不能错位。」"""
    expected = [f"C{index + 1:03d}" for index in TRUE_ORDER]
    assert solved["completed"].returncode == 0
    assert [row["coupon_id"] for row in solved["result"]["layer_order"]] == expected
    assert [row["slot_index"] for row in solved["result"]["layer_order"]] == [1, 2, 3, 4, 5]


def test_places_nonlinear_curves_in_recovered_order(solved):
    """题面：「`conductivity_curves` 与装配顺序一致，每项键集合恰好为：」"""
    expected = [PARAMETERS[index][2] for index in TRUE_ORDER]
    actual = [row["k_ref_w_mk"] for row in solved["result"]["conductivity_curves"]]
    assert solved["completed"].returncode == 0
    assert max(abs(left - right) for left, right in zip(actual, expected)) <= 1.0e-9


def test_recovers_temperature_coefficients_in_order(solved):
    """题面：「四项材料参数不得完全相同，否则对应关系不具备唯一性。」"""
    expected = [PARAMETERS[index][3] for index in TRUE_ORDER]
    actual = [
        row["temperature_coefficient_per_k"]
        for row in solved["result"]["conductivity_curves"]
    ]
    assert solved["completed"].returncode == 0
    assert max(abs(left - right) for left, right in zip(actual, expected)) <= 1.0e-9


def test_recovers_each_sensor_gain(solved):
    """题面：「传感器增益的允许偏差按 `2.0e-4` 控制。」"""
    actual = [row["gain"] for row in solved["result"]["sensor_calibration"]]
    assert solved["completed"].returncode == 0
    assert max(abs(left - right) for left, right in zip(actual, TRUE_GAINS)) <= 2.0e-4


def test_recovers_each_sensor_offset(solved):
    """题面：「传感器偏移的允许偏差按 `2.0e-4` 控制。」"""
    actual = [row["offset_k"] for row in solved["result"]["sensor_calibration"]]
    assert solved["completed"].returncode == 0
    assert max(abs(left - right) for left, right in zip(actual, TRUE_OFFSETS)) <= 2.0e-4


def test_predicts_future_sensor_temperatures(solved):
    """题面：「停热后的温度预测允许偏离实测轨迹 `0.010 K` 以内。」"""
    expected = _expected_forecast(solved["payload"])
    actual = [
        [sample["temperature_rise_k"] for sample in row["samples"]]
        for row in solved["result"]["forecast"]
    ]
    assert solved["completed"].returncode == 0
    assert max(
        abs(actual[sensor][index] - expected[index][sensor])
        for sensor in range(5)
        for index in range(len(expected))
    ) <= 0.010


def test_preserves_forecast_time_order(solved):
    """题面：「预测样行键集合恰好为：」"""
    expected = solved["payload"]["forecast"]["times_s"]
    actual = [[sample["time_s"] for sample in row["samples"]] for row in solved["result"]["forecast"]]
    assert solved["completed"].returncode == 0
    assert actual == [expected] * 5


def test_emits_exact_result_schema(solved):
    """题面：「结果根对象的键集合恰好为：」"""
    assert solved["completed"].returncode == 0
    assert set(solved["result"]) == {
        "schema_version",
        "case_id",
        "layer_order",
        "conductivity_curves",
        "sensor_calibration",
        "forecast",
        "fit",
    }


def test_emits_exact_layer_schema(solved):
    """题面：「`layer_order` 按 `slot_index` 从 `1` 到 `5` 排列。每项记录恰好包含两个字段：」"""
    assert solved["completed"].returncode == 0
    assert all(set(row) == {"slot_index", "coupon_id"} for row in solved["result"]["layer_order"])


def test_emits_exact_calibration_schema(solved):
    """题面：「`sensor_calibration` 按 `sensor_index` 从 `1` 到 `5` 排列」"""
    expected = [1, 2, 3, 4, 5]
    assert solved["completed"].returncode == 0
    assert [row["sensor_index"] for row in solved["result"]["sensor_calibration"]] == expected
    assert all(set(row) == {"sensor_index", "gain", "offset_k"} for row in solved["result"]["sensor_calibration"])


def test_emits_exact_forecast_schema(solved):
    """题面：「`forecast` 按 `sensor_index` 从 `1` 到 `5` 排列。每项记录恰好包含两个字段：」"""
    assert solved["completed"].returncode == 0
    assert all(set(row) == {"sensor_index", "samples"} for row in solved["result"]["forecast"])
    assert all(
        set(sample) == {"time_s", "temperature_rise_k"}
        for row in solved["result"]["forecast"]
        for sample in row["samples"]
    )


def test_emits_exact_fit_schema(solved):
    """题面：「`fit` 的键集合恰好为：」"""
    assert solved["completed"].returncode == 0
    assert set(solved["result"]["fit"]) == {
        "max_abs_temperature_error_k",
        "total_heat_j_m2",
        "heat_balance_error_j_m2",
    }


def test_fit_error_is_within_limit(solved):
    """题面：「训练观测的最大绝对温度误差上界为 `0.020 K`」"""
    assert solved["completed"].returncode == 0
    assert solved["result"]["fit"]["max_abs_temperature_error_k"] <= 0.020


def test_total_heat_matches_flux_integral(solved):
    """题面：「`total_heat_j_m2` 是整段输入热流对时间的积分」"""
    expected = solved["payload"]["aggregate"]["total_heat_j_m2"]
    assert solved["completed"].returncode == 0
    assert abs(solved["result"]["fit"]["total_heat_j_m2"] - expected) <= 1.0e-5


def test_heat_balance_error_is_within_limit(solved):
    """题面：「聚合热量与分段积分相差不得超过 `0.5 J/m^2`」"""
    assert solved["completed"].returncode == 0
    assert abs(solved["result"]["fit"]["heat_balance_error_j_m2"]) <= 0.5


def test_ledger_contains_both_digests(solved):
    """题面：「行对象的键集合恰好为：」"""
    record = json.loads(solved["ledger"].read_text(encoding="utf-8"))
    assert solved["completed"].returncode == 0
    assert set(record) == {
        "case_id",
        "event",
        "input_sha256",
        "result_sha256",
        "schema_version",
    }
    assert record["input_sha256"] == hashlib.sha256(solved["input_path"].read_bytes()).hexdigest()
    assert record["result_sha256"] == hashlib.sha256(solved["raw"]).hexdigest()


def test_repeated_identical_batch_does_not_append(tmp_path):
    """题面：「同一炉次拿同一热史重跑，只承认首条完成页，交班账本不能因此多出一行。」"""
    payload = _payload()
    first, _, ledger_path, _ = _invoke(tmp_path, payload)
    before = ledger_path.read_bytes()
    second, _, _, _ = _invoke(tmp_path, payload)
    assert first.returncode == 0
    assert second.returncode == 0
    assert ledger_path.read_bytes() == before
    assert len(before.splitlines()) == 1


def test_second_batch_appends_another_record(tmp_path):
    """题面：「换炉次时只在账本尾部续页。」"""
    first, _, ledger_path, _ = _invoke(tmp_path, _payload("wall-first"))
    second, _, _, _ = _invoke(tmp_path, _payload("wall-second"))
    assert first.returncode == 0
    assert second.returncode == 0
    assert len(ledger_path.read_text(encoding="utf-8").splitlines()) == 2


def test_accepts_five_layer_limit(tmp_path):
    """题面：「取值恰好等于上限时仍应处理」"""
    completed, output_path, ledger_path, _ = _invoke(tmp_path, _payload())
    assert completed.returncode == 0
    assert output_path.exists()
    assert ledger_path.exists()


def test_rejects_six_layer_rows(tmp_path):
    """题面：「任一数量或时间超过上限时拒绝」"""
    payload = copy.deepcopy(_payload())
    extra = copy.deepcopy(payload["coupons"][0])
    extra["coupon_id"] = "C006"
    extra["thickness_m"] = 0.029
    payload["coupons"].append(extra)
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_accepts_three_segment_limit(tmp_path):
    """题面：「第一段从 0.0 秒开始。」"""
    completed, output_path, ledger_path, _ = _invoke(tmp_path, _payload())
    assert completed.returncode == 0
    assert output_path.exists()
    assert ledger_path.exists()


def test_rejects_four_segments(tmp_path):
    """题面：「后一段的起点必须晚于前一段终点」"""
    payload = copy.deepcopy(_payload())
    extra = copy.deepcopy(payload["segments"][-1])
    extra["segment_index"] = 4
    extra["start_time_s"] = 16.0
    extra["end_time_s"] = 17.0
    extra["samples"][0]["time_s"] = 16.0
    extra["samples"][-1]["time_s"] = 17.0
    payload["segments"].append(extra)
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_accepts_120_samples_per_segment(tmp_path):
    """题面：「每段的 `samples` 有 48 到 120 行」"""
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path, _payload(sample_counts=(120, 120, 120))
    )
    assert completed.returncode == 0
    assert output_path.exists()
    assert ledger_path.exists()


def test_rejects_121_samples_per_segment(tmp_path):
    """题面：「每段的 `samples` 有 48 到 120 行」"""
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path, _payload(sample_counts=(121, 90, 110))
    )
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_accepts_40_forecast_points(tmp_path):
    """题面：「`times_s` 有 1 到 40 个严格递增的有限数」"""
    times = [15.3 + 0.23 * index for index in range(40)]
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path, _payload(forecast_times=times)
    )
    assert completed.returncode == 0
    assert output_path.exists()
    assert ledger_path.exists()


def test_rejects_41_forecast_points(tmp_path):
    """题面：「`times_s` 有 1 到 40 个严格递增的有限数」"""
    times = [15.3 + 0.20 * index for index in range(41)]
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path, _payload(forecast_times=times)
    )
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_top_level_extra_key(tmp_path):
    """题面：「输入根对象的键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["operator"] = "night"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_missing_segments_key(tmp_path):
    """题面：「`segments` 恰好有三个对象」"""
    payload = copy.deepcopy(_payload())
    payload.pop("segments")
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_duplicate_json_key(tmp_path):
    """题面：「读入、契约核对或反算闭合只要有一关失守」"""
    raw = b'{"schema_version":1,"schema_version":1,"case_id":"dup"}'
    completed, output_path, ledger_path = _raw_invoke(tmp_path, raw)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_boolean_thickness(tmp_path):
    """题面：「四类数值都必须是有限数，布尔值不算数值」"""
    payload = copy.deepcopy(_payload())
    payload["coupons"][0]["thickness_m"] = True
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_unsorted_forecast_times(tmp_path):
    """题面：「所有时刻都晚于最后一段终点」"""
    payload = copy.deepcopy(_payload())
    payload["forecast"]["times_s"] = [18.0, 17.0]
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_nonuniform_segment_time_axis(tmp_path):
    """题面：「段内严格等步长递增」"""
    payload = copy.deepcopy(_payload())
    payload["segments"][0]["samples"][8]["time_s"] += 0.01
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_wrong_sensor_count(tmp_path):
    """题面：「`sensor_temperatures_k` 恰好有五个有限数值」"""
    payload = copy.deepcopy(_payload())
    payload["segments"][0]["samples"][4]["sensor_temperatures_k"].pop()
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_ambiguous_material_signature(tmp_path):
    """题面：「两组卡片的四项材料参数不得完全相同」"""
    payload = copy.deepcopy(_payload())
    payload["coupons"][1]["thickness_m"] = payload["coupons"][0]["thickness_m"]
    payload["coupons"][1]["areal_heat_capacity_j_m2k"] = payload["coupons"][0]["areal_heat_capacity_j_m2k"]
    payload["coupons"][1]["k_ref_w_mk"] = payload["coupons"][0]["k_ref_w_mk"]
    payload["coupons"][1]["temperature_coefficient_per_k"] = payload["coupons"][0]["temperature_coefficient_per_k"]
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_inconsistent_total_heat(tmp_path):
    """题面：「聚合热量与分段积分相差不得超过 `0.5 J/m^2`」"""
    payload = copy.deepcopy(_payload())
    payload["aggregate"]["total_heat_j_m2"] += 10.0
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_invalid_batch_preserves_previous_files(tmp_path):
    """题面：「这个收口只留退出码，标准输出与标准错误保持零字节，已装订的层序单和交班账本也不动。」"""
    payload = copy.deepcopy(_payload())
    payload["coupons"][0]["thickness_m"] = 0.001
    old_ledger = (
        b'{"case_id":"old","event":"completed","input_sha256":"'
        + b"0" * 64
        + b'","result_sha256":"'
        + b"1" * 64
        + b'","schema_version":1}\n'
    )
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path,
        payload,
        existing_output=b"previous-layout\n",
        existing_ledger=old_ledger,
    )
    assert completed.returncode == 3
    assert output_path.read_bytes() == b"previous-layout\n"
    assert ledger_path.read_bytes() == old_ledger


def test_ledger_directory_failure_preserves_output(tmp_path):
    """题面：「任意一次换名未完成，上一版层序报告和旧账本仍须保持原字节。」"""
    payload = _payload()
    input_path = tmp_path / "case.json"
    output_path = tmp_path / "result.json"
    ledger_path = tmp_path / "ledger-dir"
    input_path.write_text(json.dumps(payload), encoding="utf-8")
    output_path.write_bytes(b"previous\n")
    ledger_path.mkdir()
    completed = subprocess.run(
        [
            "python3",
            ENTRY,
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ledger",
            str(ledger_path),
        ],
        capture_output=True,
        text=True,
        timeout=45,
        env=os.environ.copy(),
    )
    assert completed.returncode == 3
    assert output_path.read_bytes() == b"previous\n"
    assert ledger_path.is_dir()


def test_malformed_ledger_preserves_previous_output(tmp_path):
    """题面：「已归档旧页不允许回写」"""
    completed, output_path, ledger_path, _ = _invoke(
        tmp_path,
        _payload(),
        existing_output=b"older\n",
        existing_ledger=b"{broken",
    )
    assert completed.returncode == 3
    assert output_path.read_bytes() == b"older\n"
    assert ledger_path.read_bytes() == b"{broken"


def test_result_serialization_is_compact_sorted_single_line(solved):
    """题面：「层序单保持一行紧凑布局，所有键名按码点递增排列。」"""
    parsed = json.loads(solved["raw"])
    assert solved["completed"].returncode == 0
    assert b", " not in solved["raw"]
    assert solved["raw"].endswith(b"\n")
    assert not solved["raw"].endswith(b"\n\n")
    assert list(parsed) == sorted(parsed)


def test_success_path_is_silent(solved):
    """题面：「成功时标准输出、标准错误两个流都保持零字节。」"""
    assert solved["completed"].returncode == 0
    assert solved["completed"].stdout == ""
    assert solved["completed"].stderr == ""


def test_case_identifier_is_preserved(solved):
    """题面：「`case_id` 是非空字符串」"""
    assert solved["completed"].returncode == 0
    assert solved["result"]["case_id"] == "wall-nonlinear-0611"


def test_limits_object_rejects_extra_key(tmp_path):
    """题面：「`limits` 的键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["limits"]["units"] = "si"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_input_declares_five_calibrated_channels(tmp_path):
    """题面：「五个热电偶各自换过放大器」"""
    payload = _payload()
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 0
    assert output_path.exists()
    assert all(
        len(sample["sensor_temperatures_k"]) == 5
        for segment in payload["segments"]
        for sample in segment["samples"]
    )
    assert ledger_path.exists()


def test_rejects_wrong_schema_version(tmp_path):
    """题面：「`schema_version` 是整数 `1`。」"""
    payload = copy.deepcopy(_payload())
    payload["schema_version"] = 2
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_extra_coupon_key(tmp_path):
    """题面：「`coupons` 恰好有五个对象。每个对象的键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["coupons"][0]["vendor"] = "north"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_extra_sample_key(tmp_path):
    """题面：「每段的 `samples` 有 48 到 120 行。样行键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["segments"][0]["samples"][0]["unit"] = "K"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_extra_aggregate_key(tmp_path):
    """题面：「`aggregate` 的键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["aggregate"]["unit"] = "J/m^2"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_wrong_sample_count(tmp_path):
    """题面：「`sample_count` 是三段有效样本行数之和，且必须与 `samples` 实际行数一致。」"""
    payload = copy.deepcopy(_payload())
    payload["aggregate"]["sample_count"] += 1
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_rejects_extra_forecast_key(tmp_path):
    """题面：「`forecast` 的键集合恰好为：」"""
    payload = copy.deepcopy(_payload())
    payload["forecast"]["unit"] = "s"
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_calibration_rows_have_three_exact_fields(tmp_path):
    """题面：「`sensor_calibration` 按 `sensor_index` 从 `1` 到 `5` 排列。每项记录恰好包含三个字段：」"""
    completed, output_path, ledger_path, _ = _invoke(tmp_path, _payload())
    result = json.loads(output_path.read_bytes())
    assert completed.returncode == 0
    assert all(
        set(row) == {"sensor_index", "gain", "offset_k"}
        for row in result["sensor_calibration"]
    )
    assert ledger_path.exists()


def test_rejects_segment_end_after_thirty_seconds(tmp_path):
    """题面：「记录窗口的公开上限是 30.0 秒。」"""
    payload = copy.deepcopy(_payload())
    payload["segments"][-1]["end_time_s"] = 30.5
    payload["segments"][-1]["samples"][-1]["time_s"] = 30.5
    completed, output_path, ledger_path, _ = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()
    assert not ledger_path.exists()


def test_ledger_rows_are_compact_single_lines(tmp_path):
    """题面：「每条完成记录保持一行紧凑布局，键名按码点递增排列，记录内部不留格式空白，末尾带一个换行字节。」"""
    completed, output_path, ledger_path, _ = _invoke(tmp_path, _payload())
    output_raw = output_path.read_bytes()
    raw = ledger_path.read_bytes()
    assert completed.returncode == 0
    assert b", " not in output_raw
    assert output_raw.endswith(b"\n")
    assert not output_raw.endswith(b"\n\n")
    assert b", " not in raw
    assert raw.endswith(b"\n")
    assert not raw.endswith(b"\n\n")
    assert len(raw.splitlines()) == 1


def test_full_batch_completes_within_public_limit(tmp_path):
    """题面：「单批次从读入到写出应在 40 秒内完成。」"""
    payload = _payload()
    input_path = tmp_path / "timed-input.json"
    output_path = tmp_path / "timed-result.json"
    ledger_path = tmp_path / "timed-ledger.jsonl"
    input_path.write_text(json.dumps(payload), encoding="utf-8")
    started = time.monotonic()
    completed = subprocess.run(
        [
            "python3",
            ENTRY,
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ledger",
            str(ledger_path),
        ],
        capture_output=True,
        text=True,
        timeout=45,
        env=os.environ.copy(),
    )
    elapsed = time.monotonic() - started
    assert completed.returncode == 0
    assert elapsed < 40.0
