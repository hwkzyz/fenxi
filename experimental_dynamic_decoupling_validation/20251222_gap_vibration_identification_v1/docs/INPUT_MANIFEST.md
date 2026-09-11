# 随附输入清单

本版本随附下列正式输入，因而 Main09/Main10 不需要读取旧程序目录：

- `OPR6_SlotAngleCalibration_Laps2_21_20251222.mat/csv`
- `Step05_Response_Surface_20251222.mat`
- `Step05I_OffsetTilt_Shared_Response_Surface_20251222.mat`
- `Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat`
- `Step06I_OffsetTiltShared_GapLibrary_20251222_B5_S123.mat`
- `inputs/prepared/foundation/step01_low_speed_reference`
- `inputs/prepared/foundation/step02_dynamic_btt`
- `inputs/prepared/foundation/step03_observation_bundle`
- `inputs/prepared/foundation/step04_low_speed_template`

原始传感器和应变数据仍是外部实验数据，由 `Config_20251222.m` 统一指向
`E:\试验数据\20251222`。静态间隙原始库和叶片间隙统计表也是数据源，只有
重新运行 Main05/Main06 时才需要。
