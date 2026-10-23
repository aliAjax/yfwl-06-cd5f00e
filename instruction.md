# 多频热脉冲下的五层壁复核

夜班复测留下的五层炉壁记录准备交给材料组。层卡没有按装配顺序归档，五个热电偶各自换过放大器，三段不同频段的热脉冲之间还隔着短暂停热。现在需要从观测温度史和整段累计热量反推装配关系、传感器刻度，并给出停热后的短时温度预测。

入口固定为 `/app/recover_thermal.py`，调用方式：

```text
python3 /app/recover_thermal.py --input <输入路径> --output <结果路径> --ledger <账本路径>
```

单批次从读入到写出应在 40 秒内完成。

输入根对象的键集合恰好为：

```text
schema_version, case_id, coupons, segments, aggregate, forecast, limits
```

`schema_version` 是整数 `1`。`case_id` 是非空字符串。

`coupons` 恰好有五个对象。每个对象的键集合恰好为：

```text
coupon_id, thickness_m, areal_heat_capacity_j_m2k, k_ref_w_mk, temperature_coefficient_per_k
```

五个 `coupon_id` 互不相同。`thickness_m` 位于 `0.008` 到 `0.080` 米之间。`areal_heat_capacity_j_m2k` 位于 `400.0` 到 `5000.0` 之间。`k_ref_w_mk` 位于 `0.20` 到 `10.0` 之间。`temperature_coefficient_per_k` 位于 `-0.050` 到 `0.050` 之间。四类数值都必须是有限数，布尔值不算数值。两组卡片的四项材料参数不得完全相同，否则对应关系不具备唯一性。

`segments` 恰好有三个对象，键集合恰好为：

```text
segment_index, frequency_hz, phase_rad, start_time_s, end_time_s, samples
```

`segment_index` 从 `1` 开始连续编号。`frequency_hz` 位于 `0.03` 到 `0.50` 之间。`phase_rad` 位于 `0.0` 到 `6.283185307179586` 之间。第一段从 0.0 秒开始。后一段的起点必须晚于前一段终点。记录窗口的公开上限是 30.0 秒。

每段的 `samples` 有 48 到 120 行。样行键集合恰好为：

```text
time_s, heat_flux_w_m2, sensor_temperatures_k
```

`time_s` 从该段 `start_time_s` 开始，到该段 `end_time_s` 结束，段内严格等步长递增。`heat_flux_w_m2` 位于 `1.0` 到 `2000.0` 之间。`sensor_temperatures_k` 恰好有五个有限数值，依次属于界面传感器 `1` 到 `5`。段外热流为零；段内热流按相邻样行做线性插值。

`aggregate` 的键集合恰好为：

```text
total_heat_j_m2, sample_count
```

`total_heat_j_m2` 是整段输入热流对时间的积分，单位 `J/m^2`。`sample_count` 是三段有效样本行数之和，且必须与 `samples` 实际行数一致。聚合热量与分段积分相差不得超过 `0.5 J/m^2`。

`forecast` 的键集合恰好为：

```text
times_s
```

`times_s` 有 1 到 40 个严格递增的有限数，所有时刻都晚于最后一段终点，最远不超过最后一段终点后 `10.0` 秒。预测阶段热流为零。

`limits` 的键集合恰好为：

```text
max_layers, max_coupons, max_segments, max_samples_per_segment, max_forecast_points, max_time_s
```

各项依次固定为 `5`、`5`、`3`、`120`、`40` 和 `30.0`。取值恰好等于上限时仍应处理；任一数量或时间超过上限时拒绝。

第 `i` 个界面节点的温升 `T_i(t)` 按五节点热网络演化。令 `T_0=0`，第 `i` 层的导热系数为：

```text
k_i(T) = k_ref_i * (1 + temperature_coefficient_i * T)
g_i(T) = k_i(T) / thickness_i
```

每层的界面热容取该层的 `areal_heat_capacity_j_m2k`。相邻节点的热导按两侧节点温升的平均值计算。节点 `1` 到 `4` 的热平衡为：

```text
C_i dT_i/dt = g_i(T_avg_i) * (T_(i-1) - T_i) + g_(i+1)(T_avg_(i+1)) * (T_(i+1) - T_i)
```

节点 `5` 的热平衡为：

```text
C_5 dT_5/dt = g_5(T_avg_5) * (T_4 - T_5) + q(t)
```

五个传感器分别贴合五个界面。观测关系为：

```text
observed_s(t) = gain_s * T_s(t) + offset_s
```

`gain` 位于 `0.70` 到 `1.30` 之间，`offset_k` 位于 `-0.30` 到 `0.30` 之间。训练观测的最大绝对温度误差上界为 `0.020 K`。

结果根对象的键集合恰好为：

```text
schema_version, case_id, layer_order, conductivity_curves, sensor_calibration, forecast, fit
```

`layer_order` 按 `slot_index` 从 `1` 到 `5` 排列。每项记录恰好包含两个字段：

```text
slot_index, coupon_id
```

`conductivity_curves` 与装配顺序一致，每项键集合恰好为：

```text
slot_index, k_ref_w_mk, temperature_coefficient_per_k
```

`sensor_calibration` 按 `sensor_index` 从 `1` 到 `5` 排列。每项记录恰好包含三个字段：

```text
sensor_index, gain, offset_k
```

`forecast` 按 `sensor_index` 从 `1` 到 `5` 排列。每项记录恰好包含两个字段：

```text
sensor_index, samples
```

预测样行键集合恰好为：

```text
time_s, temperature_rise_k
```

预测值要按本路传感器刻度回算。装配序列不能错位。传感器增益的允许偏差按 `2.0e-4` 控制。传感器偏移的允许偏差按 `2.0e-4` 控制。停热后的温度预测允许偏离实测轨迹 `0.010 K` 以内。

`fit` 的键集合恰好为：

```text
max_abs_temperature_error_k, total_heat_j_m2, heat_balance_error_j_m2
```

`max_abs_temperature_error_k` 是训练观测上的最大绝对拟合误差。`total_heat_j_m2` 是按输入热流重新积分得到的累计热量。`heat_balance_error_j_m2` 是输入聚合热量减去该积分值。

层序单保持一行紧凑布局，所有键名按码点递增排列。逗号和冒号之间不插空白，文件结尾保留一个换行字节。热工交接要求控制台无提示音也无文字，成功时标准输出、标准错误两个流都保持零字节。

每次成功都要在 `--ledger` 追加一条记录。行对象的键集合恰好为：

```text
case_id, event, input_sha256, result_sha256, schema_version
```

`event` 是 `completed`。输入文件和结果文件分别计算小写 SHA-256。同一炉次拿同一热史重跑，只承认首条完成页，交班账本不能因此多出一行。已归档旧页不允许回写，换炉次时只在账本尾部续页。每条完成记录保持一行紧凑布局，键名按码点递增排列，记录内部不留格式空白，末尾带一个换行字节。

读入、契约核对或反算闭合只要有一关失守，程序就按 `3` 号结束信号收口。这个收口只留退出码，标准输出与标准错误保持零字节，已装订的层序单和交班账本也不动。结果与账本先各自落到同目录的炉次暂存件，再依次换名。任意一次换名未完成，上一版层序报告和旧账本仍须保持原字节。
