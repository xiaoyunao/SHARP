# 隔离的新模式关联原型

仅供本地MP晨昏模式。复制基线为954e43d5的make_tracklet_linreproj.py，并保留380b777b验证过的实际曝光中点/HDU时间适配；不导入、不修改原heliolincrr。时间helper与关联器均位于本目录。直接运行须显式给出本地root/outdir，无服务器默认路径。

推荐入口：`python3 -m targeted_observing.offline_processing --night YYYYMMDD --root /local/data --output /local/output`，速度范围由新模式配置给出3–750″/小时。原生产代码保持原样。此原型没有MPC提交功能，也不代表完整unknown管线已端到端验证。

后续改动只在本目录进行；生产代码若独立演进，需要显式评审是否迁移，不能直接覆盖旧目录。
