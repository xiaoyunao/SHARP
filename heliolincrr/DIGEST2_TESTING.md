# digest2 独立测试工具

这是接入前的独立测试实现；daily、exporter、review packager、watcher 和 survey
均没有接入新模块。不能用这套测试账本代替正式 MPC submission 状态。

## 安装与运行（本机或独立临时环境）

```bash
python3 -m venv /tmp/smt-digest2-test-env
/tmp/smt-digest2-test-env/bin/python -m pip install -r heliolincrr/requirements-digest2-test.txt
/tmp/smt-digest2-test-env/bin/python -m unittest discover -s heliolincrr/tests -p 'test_digest2_rehearsal.py' -v
```

程序要求 Python >=3.10，固定 `digest2==2.8.2`。该版本发布的预编译 wheel
只有 CPython 3.12；其他 Python 版本可能需要从源码编译。不要在生产环境内
直接运行安装命令。模型及 MPC.config 随 digest2 包提供；ObsCodes 需另行缓存。

评分输入必须是从生产数据只读复制的 PSV、与该 PSV 对象集合匹配的 catalog，
以及可选的 review CSV。回放历史上报用 `submit_masked.json`；回放全部未知源
用完整 unknown catalog 和未经人工过滤的观测导出。

```bash
/tmp/smt-digest2-test-env/bin/python heliolincrr/score_unknown_digest2.py \
  --psv /path/to/copied_observations.psv \
  --catalog /path/to/matching_unknown_catalog.json \
  --review-csv /path/to/copied_review.csv \
  --obscodes /path/to/ObsCodes.html \
  --out-dir /path/to/new_rehearsal_directory
```

- 输出目录必须尚不存在，防止覆盖之前的输入/结果。
- 不传 review CSV 时仍为所有输入目标评分，但不会生成待验证批次。
- 保留原始 ADES 时间/坐标/误差/星等；不会自动给历史时间加 15 秒。
- 采用官方 PSV parser，保留 G -> V 的星等转换，使用 ADES 误差模式和默认
  MPC 误差下限/上限。2.8.2 已支持 PSV，且程序额外验证行数，防止其 parser
  静默跳过坏行后对残缺轨迹评分。
- `routing.json` 保存逐目标分数、真假决定、分流、错误、输入和依赖 hash。
  `inputs/` 保存本次评分快照；模型、MPC.config、ObsCodes 也被快照。
- 人工接受且 NoID NEO >=65 进入 `neo_candidate.psv`；低分进入
  `unclassified.psv`；假源、未复核和评分失败对象不进入批次。
- catalog 缺失观测、重复历元、错误台址、非有限值和已编号对象被拒绝。
  当前仅支持单台址、单夜、地基 CCD、项目 G 波段和单个 ADES block。
- CLI 在存在评分失败时返回 2，同时保存可诊断结果；不能忽略退出码。

## 验证与状态测试

默认只做本地完整性检查，没有网络请求：

```bash
/tmp/smt-digest2-test-env/bin/python heliolincrr/validate_neo_batches.py \
  /path/to/new_rehearsal_directory/routing.json
```

只有显式使用 `--test-endpoint --email <address>` 才调用 MPC
`https://minorplanetcenter.net/submit_psv_test`。程序不提供正式提交 URL 或
`--submit` 参数。输入是内存中已校验 hash 的字节，通过 `source=<-` 发送。

逐批状态保存在结果目录的 `validation_state.json`；单元测试的模拟状态保存
在独立 `simulation_state.json`，两个命名空间不能混用。

- 跨进程锁覆盖重新加载状态、检查、请求和原子/fsync 落盘。
- 每批验证成功即落盘，下一批失败不影响成功批；重启不重发已成功批。
- HTTP 200 / CurlID 不视为格式验证通过；须有明确的格式成功回复。
- 只有 Submission ID 的回执记录为 `test_received`、`format_confirmed=false`，
  不判为拒绝、不自动重发，也不判为格式通过；CLI 仍返回 2，避免调用方只凭退出码放行。
- 明确 HTTP 4xx 记录为拒绝；超时、5xx 或不认识的回复保留为待核对。
- 超时/异常/中断的请求标记 `needs_reconciliation`，不会盲目自动重试。
- 被 MPC 明确拒绝的批次保存 `rejected`，修复后生成新快照/新批次。
- 在第一个请求前检查所有批次：输入及 payload hash、评分与真假门控、
  观测集合、分数阈值、批次不重叠且并集完整。

这些状态仅验证测试工具行为。上线时仍需实现正式逐批提交账本、历史状态迁移、
NEOCP/已知身份路由，并确认正式回执解析；不能直接把 `validated` 当作 `submitted`。

## 2026-09-20 实际回放结果

使用当时保存的上报 PSV，保持历史时间约定，NoID NEO 分数如下。
这是使用冻结模型/配置的本地重算，不代表 MPC 内部评分的精确副本，也不确认目标真实性。

| trkSub | 夜次 | 观测数 | NoID NEO | Raw NEO |
| --- | --- | ---: | ---: | ---: |
| 00001nF | 20260902 | 3 | 70.9705 | 40.6510 |
| 00001nP | 20260902 | 3 | 99.4856 | 90.5651 |
| 00001o2 | 20260902 | 3 | 99.7659 | 90.7988 |
| 00001oI | 20260902 | 3 | 99.8908 | 96.3044 |
| 00001oj | 20260902 | 4 | 73.5258 | 16.9152 |
| 00001ua | 20260910 | 3 | 99.9713 | 99.6494 |
| 00001vy | 20260910 | 3 | 99.2807 | 67.7763 |

- 两夜实际上报共 9 个目标（另含 `00001tV`, `00001tY`），与官方
  `Digest2.classify_file(..., max_workers=1)` 的 Raw/NoID 结果逐项完全一致。
- 20260910 全部 120 条 unknown / 372 条观测评分成功：107 条 >=65、13 条 <65。
  三次不同复核输入的运行得到完全相同的 120 个 NoID 分数。
- 使用真实历史 submit CSV：只有 4 条进入 NEO 候选批，116 条被人工门控排除。
- 为覆盖两个非空批次，另用明确标记 `SIMULATED_all_accept` 的假定复核 CSV
  做拆批压力测试：107/13 两批，集合不重叠、观测总数守恒；这不是真实复核结论。
- 按 09-04 WORKLOG 的四个假源决定做单独门控测试：四条均排除，只有
  `00001nF` 保留；这不构成对 `00001nF` 的新增真实性确认。
- 官方 MPC80 教程样例成功解析评分；该样例与教程中 programmatic 样例的
  坐标不同，不混用两套数据的预期结果。
- 21 项单元/并发/CLI 测试通过；测试传输均为 mock，没有向 MPC POST。

服务器基线：Python 3.10.20 / x86_64 / glibc 2.39，有 gcc；原环境未安装
digest2。生产 exporter/watch_submit_reviews 的 SHA-256 在只读回放前后相同。

服务器隔离目录：`/tmp/smt-digest2-rehearsal-20260920.4ihjJN`。由于服务器直连
PyPI 下载很慢，改为本机下载 sdist/依赖 wheel、通过 scp 送入该目录后离线构建。
安装成功于独立 venv；生产 heliolinc 的 `find_spec('digest2')` 仍为 `None`。
同一套 21 项测试也在服务器 Python 3.10.20 通过；两夜共 125 个目标的
Linux/macOS 分流及模型/配置/ObsCodes hash 一致，分数差异小于 `1e-8`。

本地完整结果目录：`outputs/digest2_rehearsal_20260920/`（不提交 git）。
包含 `mpc_seven_scores.csv`、`all_20260910_scores.csv`、输入 hash 清单、
各场景 `routing.json`/PSV/依赖快照、双平台测试日志和环境版本，以及源码安装包。

## 2026-09-20 MPC 测试接口实际回执

经用户授权，从本机向 `https://minorplanetcenter.net/submit_psv_test` 发出
4 次测试 POST，回执地址使用现有配置 `wsgp2024@163.com`。正式端点请求数为 0。
服务器生产代码、数据、环境、定时任务均未修改。

两个正常样本分别只有 3 条观测；低分样本 `00001u6` 的真实历史复核为否，
本次仅用明确标记 FORMAT_TEST_ONLY 的假定复核 CSV 测试 ordinary 分支格式，
没有更改其真实复核结论。测试结果：

| 测试 | 即时回复 |
| --- | --- |
| `neo candidate`，00001ua | `2026-09-20T03:52:14.232_00000kOY` |
| `unclassified`，00001u6 | `2026-09-20T03:52:16.547_00000kOZ` |
| 故意过短的无效文本 | HTTP 400：data size 不在 200..1000000000 bytes 范围 |
| 保留合法 header/长度，只把首行 RA 改为 NOT_A_NUMBER | `2026-09-20T03:52:59.895_00000kOa` |

三份 ID 回执正文都只有 ACK 文本和 `Submission ID is ...`，没有
`Submission format valid.`。由于错误 RA 也能收到 ID，不能据此声称 ADES
字段级格式验证通过。官方命令行文档也说明测试接口返回 CurlID，且不会把测试批次
送入正式处理：
https://docs.minorplanetcenter.net/mpc-ops-docs/observations/command-line-submissions/

已修正测试回执分类器：收件编号单独记为 `test_received`；对现有两条记录仅离线
重判，没有重发。此前误标 rejected 的原账本已保留。加入真实回执与 400 对照的
回归测试，本机 25/25 PASS。原始回执、请求样本和账本见
`outputs/digest2_rehearsal_20260920/mpc_format_test/`。

剩余上线工作：统一未来观测的曝光中点时间、复核页面分数展示、把评分门控接到
所有正式入口、逐批正式回执和历史去重迁移；另外建立可证明的 ADES schema/字段
校验，不能以测试 CurlID 作为通过证据。
本次不修改任何生产入口，也不重新上报历史观测。
