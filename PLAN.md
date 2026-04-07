# PLAN

## Current objective

把 `smt_asteroid` 整理成清晰的双版本仓库：

- `*_server` 保存服务器现版备份
- `*_local` 保存本地可运行版本

## Milestones

1. 完成目录重组
2. 补齐项目总 README 和项目级 `AGENTS.md`
3. 清理旧版本目录和运行产物
4. 做一轮最小验证
5. 初始化 git 并连接 GitHub 远端

## Outstanding issues

- `known_asteroid_server` 的非核心辅助脚本尚未全部从服务器回收
- 本地大依赖文件不应进入 git
- GitHub 远端仓库名称/目标尚未最终确定

## Validation criteria

- 新目录结构只保留 `*_local` 和 `*_server`
- 本地脚本默认路径不再指向服务器绝对路径
- 服务器备份脚本保留服务器路径语义
- 核心 Python 文件通过最小语法检查
- 已完成 `git init`

## Next recommended steps

1. 检查并提交当前整理结果
2. 确认并配置 GitHub 远端
3. 如需补齐服务器辅助脚本，再次从服务器回收
