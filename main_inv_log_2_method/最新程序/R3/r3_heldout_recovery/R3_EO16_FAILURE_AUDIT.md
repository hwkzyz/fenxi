# R3 自适应法 EO16 失败复核

日期：2026-08-27。

正式矩阵的 2700 个自适应非零振动样本中有 8 个失败，均把真实 EO10 识别为 EO16：

- `g=0.9 mm, A=0.20 mm`：相位/种子配对编号 69、72、74；
- `g=0.9 mm, A=0.25 mm`：编号 4、74；
- `g=1.0 mm, A=0.20 mm`：编号 72、73、74。

`Run_12_EOCompetitionDeltaCoordinate` 使用低速锚点一致的状态坐标
`g0+(gTruth-0.8 mm)`，重新比较每个样本的 EO10 与 EO16。

结果为 8/8 `joint_search_rank_inversion_confirmed`：

- 自由间隙搜索时，EO16 相对 EO10 的 RMSE 优势为 `3.24–385.29 μV`；
- 固定到增量一致状态后，EO10 相对 EO16 的 RMSE 优势为 `514.82–876.99 μV`。

因此，这 8 个失败不是 8 种噪声强度，也不是系统真实出现 EO16。75 表示 75 个确定性相位与随机
种子的配对记录；失败集中在其中少数记录。其机理是特定相位/噪声实现下，自由间隙参数与错误阶次
共同降低残差，使 EO16 在联合搜索中暂时优于真实 EO10。固定正确状态后排序全部恢复。

可追溯文件：

- `Run_12_EOCompetitionDeltaCoordinate.m`
- `output/12_eo_competition_delta_coordinate/eo_competition_delta_detail.csv`
- `output/12_eo_competition_delta_coordinate/eo_competition_delta_summary.csv`
