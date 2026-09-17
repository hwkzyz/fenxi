# R4 仿真冻结记录

冻结日期：2026-09-05

## 冻结主链

`Main_01_BuildSimulationSurface` -> `Main_02_PC1SurfaceValidation` ->
`Main_03_OpenSetTransfer` -> `Main_04_BoundaryIdentifiability` ->
`Main_05_PC1Figures`

## 冻结结果

- PC1 explained ratio：95.99%。
- PC1 transfer mean/worst NRMSE：0.1652% / 0.3975%。
- Oracle gap-law floor：0.0871%。
- Discrete direct baseline：0.5747%。
- Anchor-preserving incremental baseline：0.5317%。
- Interpolation cases：30/42，mean NRMSE 0.1262%。
- Endpoint extrapolation cases：12/42，mean NRMSE 0.2628%，worst 0.3975%。
- Main04 known/unknown-anchor mean NRMSE：0.1717% / 0.1745%。

## 方法边界

该结果证明的是：在当前 coefficient representation 下，单条 anchor 可以定位并迁移
目标 response surface。latent 坐标和 nuisance gap 参数不被解释为唯一物理量；Main04
的联合搜索也不构成绝对机械间隙真值验证。

## 冻结约束

R4 不再新增 Main 程序、不再更换 gap-law、`gref`、PCA/normalization 或动态求解器。
后续工作转入 R5，逐节点核查实验侧是否实现相同的 `q -> PC1 -> q(z) -> anchor/profile`
前端，并将结果送入冻结的 Main10。
