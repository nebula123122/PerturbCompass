# PerturbCompass

本项目研究：在目标细胞背景只测量少量 anchor 扰动时，如何选择或平均透明的迁移校准模型，并评估这种决策是否可靠。

当前论文主线是 APMA：对 anchor 做交叉验证，按 bootstrap 中各家族胜出的频率决定单选或概率平均。主评价采用相对最佳固定候选的 regret；另外检验重复测量一致性、source 功能簇外推与 anchor 预算。

**目录与入口**

- `configs/`：按研究用途命名的实际运行协议
- `src/perturbcompass/`：正式方法包和论文分析入口
- `results/reference/`：锁定结果，保持原始字节和哈希
- `data/processed/`：固定分区、资格列表与 source-only 设计
- `manuscript/`：当前方法、结果和待补事项
- `docs/research_history.html`：早期探索、失败、修订及完整原始记录，支持全文搜索
- `docs/provenance/`：必要冻结证据、旧新路径对照和本次整理审计

GitHub自动检查会在Linux和Windows上运行测试、证据校验及源码打包，不下载原始实验数据。服务器正式工作目录现在直接对应已发布的`main`分支，操作见 [Git与服务器说明](docs/repository_operations.md)。

安装及使用见 [英文README](README.md) 与 [复现说明](docs/reproducibility.md)。测试和从结果表复算主结论不需要下载大体积原始数据。模型重跑写入 `results/reproduced/`，不会覆盖参考结果。

当前目录中已去除历史探索runner、大体积无关数据、下载碎片和缓存。完整恢复归档位于主项目外；其私有位置记录不进入发布包。原始冻结记录正文仍保留历史编号，以便核对当时的协议与哈希。

论文所需的负结果仍保留，包括 Adamson 主门槛未通过、Nadig 双屏选择器失败、CD4 T与重复测量的盲态可行性修订，以及 anchor 预算无效首轮。这些记录不能因目录整理而从研究叙事中消失。

尚需作者确认作者信息、代码许可证和最终论文引用。源码已上传至 [nebula123122/PerturbCompass](https://github.com/nebula123122/PerturbCompass)。原始数据与恢复归档留在服务器，不进入GitHub仓库。
