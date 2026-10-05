# 利用结构预测辅助设计 RfxCas13d collateral activity 增强突变

本项目介绍如何结合实验结构、蛋白质结构预测和突变效应评估，提出并筛选可能增强 RfxCas13d（CasRx）旁切活性（collateral / trans-cleavage activity）的候选突变。

**结构预测用于提出和筛选假设，不能直接预测旁切活性倍数，也不能证明突变增强。** 最终结论需要在相同条件下比较 WT 与突变体的靶 RNA 激活依赖性及非靶 RNA 切割动力学。

资料与工具文档查阅日期：2026-10-06。本文是使用指南；尚未执行结构预测、确定增强突变或完成实验验证。

## 1. 先定义要增强的指标

这里的 collateral activity 指：crRNA–靶 RNA 识别激活后，Cas13d 对不与该 crRNA 配对的 RNA 的切割。它与错配靶 RNA 的识别/切割是不同指标。

建议把研究目标写为：在等量有效蛋白、相同 crRNA、靶 RNA 和报告 RNA 条件下，提高靶 RNA 激活后的非靶 RNA 初始切割速率，同时保留靶向活性，且不明显增加无靶背景切割。

细胞内 collateral effects 还受靶转录本丰度、蛋白表达量等因素影响。报告信号更强或细胞损伤更大，不足以单独说明蛋白本征旁切能力增强。[人细胞 RfxCas13d 旁切研究](https://www.nature.com/articles/s42003-023-04708-2)

## 2. 本目录已有资料与序列风险

| 文件 | 用途 | 使用前检查 |
| --- | --- | --- |
| `Cas13d/P54977_全序列图谱.dna` | SnapGene 质粒图谱 | 导出实际 Cas13d CDS 的翻译序列；确认背景突变及融合标签 |
| `Cas13d/9LZT_B.prot` | SnapGene 蛋白序列及注释 | 这是二进制 `.prot`，不是 FASTA；需从 SnapGene 导出 |
| `Cas13d/readme.txt` | 已有序列来源备注 | 包含 hfCas13d/N2V8 来源与参考序列一致的描述，需重新比对确认 |
| `High-fidelity Cas13 variants for targeted RNA degradation with minimal collateral effects/` | 高保真变体论文 | 该工作的主要目标是降低旁切，不能直接当作增强结果 |
| `Cas13d_Cas12a2_反式切割活性增强调研.pdf`、`20260921.pptx` | 既有调研 | 将 RfxCas13d、其他 Cas13d 同源物及 Cas12a2 的结果分别处理 |

特别注意：**N2V8 高保真背景与 WT 背景必须分别命名和比较。** 在高保真背景上恢复旁切，不等于超过 WT 的旁切活性。原论文提供可用于研究旁切调控的线索，但不能假定把降低旁切的设计“反过来”就会增强。[高保真 Cas13 原论文](https://www.nature.com/articles/s41587-022-01419-7)

现有备注还提到 R239/H244/R858/H863 的催化失活替换，并使用 `9LZU_B` 的称呼。当前 RCSB 的 9LZT 蛋白链是 B；9LZU 蛋白链是 A，两条目摘要均标注“Mutation(s): No”。这些描述存在需要核对的差异，不能只凭文件名或网页摘要决定 WT 身份。应以论文构建信息、mmCIF 聚合物序列、实际 CDS 三者比对为准。[9LZT](https://www.rcsb.org/structure/9LZT)、[9LZU](https://www.rcsb.org/structure/9LZU)

## 3. 选择工具

| 工具 | 在本流程中解决的问题 | 不能据此得出的结论 |
| --- | --- | --- |
| ColabFold / AlphaFold2 | 预测 WT 和突变体蛋白整体折叠 | 单体结构不能充分描述 RNA 激活状态；高置信度不代表高活性 |
| AlphaFold 3 / AlphaFold Server | 建模蛋白–crRNA、蛋白–crRNA–靶 RNA 复合物 | 复合物预测不等于结合自由能、激活概率或切割速率 |
| FoldX | 在固定结构附近评估替换后的能量与局部适配 | 稳定性 ΔΔG 不等于旁切活性变化 |
| PyMOL / ChimeraX / RCSB Mol* | 比对模型、查看侧链、RNA 接触与结构可信度 | 视觉上更开放或更接近 RNA，不足以证明增强 |

先用 ColabFold 做少量蛋白模型检查，再对有明确机制假设的候选做 RNA 复合物分析。ColabFold/AF2 的蛋白多聚体模式不能用来替代蛋白–RNA 建模。工具入口：[ColabFold](https://github.com/sokrypton/ColabFold)、[AlphaFold 3](https://github.com/google-deepmind/alphafold3)、[FoldX BuildModel](https://foldxsuite.crg.eu/command/BuildModel)。

## 4. 准备序列与实验结构

### 4.1 确认工作序列

1. 从 SnapGene 导出实际 CDS 的翻译序列，单独保存 Cas13d 主体序列；同时保留完整表达构建序列。
2. 与参考 [NCBI SFX39573.1](https://www.ncbi.nlm.nih.gov/protein/SFX39573.1) 和 [Addgene #109049](https://www.addgene.org/109049/) 的构建说明核对。
3. 列出所有氨基酸差异、起始/终止位置、NLS、连接肽、EGFP 等融合片段。不要将 `.dna` 文件名中的 `P54977` 直接解释为已确认的 RfxCas13d UniProt accession。
4. 确认后保存为 `inputs/RfxCas13d_WT.fasta`。如果实际是 N2V8，就保存为 `RfxCas13d_N2V8.fasta`，另外准备真实 WT。
5. 建立 `inputs/residue_mapping.csv`，记录 `WT_position,WT_aa,construct_position,pdb_id,chain_id,auth_seq_id,label_seq_id,observed`。

FASTA 只能包含真正的氨基酸序列；示意占位符不能提交预测。缺失坐标不等于序列缺失，不能删除未解析残基后继续按原编号设计突变。

### 4.2 选择结构参照

优先使用 RfxCas13d 自身结构：9LZT 是 crRNA 二元复合物；9LZU 的关联论文描述捕获短 off-target RNA 的三元复合物，不能直接当作完全匹配靶 RNA 激活的催化状态。下载 mmCIF 和验证报告，检查链编号、缺失区段及局部密度。[9LZT](https://www.rcsb.org/structure/9LZT)、[9LZU](https://www.rcsb.org/structure/9LZU)

6E9F 是 **EsCas13d** 的 crRNA–靶 RNA 三元复合物。可作同源机制参照，但必须先做序列和结构比对，不能直接把它的残基编号移到 RfxCas13d。[6E9F](https://www.rcsb.org/structure/6E9F)

## 5. 先提出机制假设，再选择突变

结构预测工具通常接受“已指定的突变序列”，并不会自动给出增强旁切的可靠位点。候选选择应结合实验结构、保守性、文献突变结果和局部结构。

以下是研究假设，均不是已验证的增强策略：

| 假设方向 | 应查看的结构特征 | 需要排除的另一种解释 |
| --- | --- | --- |
| 改善非靶 ssRNA 接近催化表面的机会 | 表面电荷、空间阻挡、溶剂暴露及 RNA 接触 | RNA 结合过强可能妨碍周转；新增接触也可能非特异吸附 |
| 改变 RNA 识别向 HEPN 激活的构象耦联 | 域间接触、连接区、多个状态之间的差异 | 可能仅放宽错配识别，或降低匹配靶活性 |
| 提高活性蛋白比例或稳定性 | 埋藏侧链、局部堆积、聚集倾向 | 可溶表达增加与单个活性酶的催化增强需要分开测量 |

首轮建议选择少量、具有独立机制理由的单点替换，先建立单点证据，再考虑组合。保护催化基序和关键 crRNA 结构支撑位点；催化失活变体可作为阴性对照。

对每个候选填写：残基映射、文献依据、预期影响、结构证据、可能的不利影响、验证指标。其他 Cas13 亚型、同源物或专利中的位点，只能在核对背景和证据后作为线索。

## 6. 使用 ColabFold 比较 WT 与突变体

### 6.1 网页方式：先验证 WT

打开 [官方 AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)，选择 GPU runtime，输入已确认的 WT 蛋白序列，运行并下载完整结果。对突变体使用相同模型设置；保存输入、MSA、置信度文件和模型，避免只保存截图。

Windows 本地批处理可使用 WSL2 下的 [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)；安装按其当前说明执行。近千残基蛋白的显存需求取决于设置，先运行一个 WT 作资源检查。

### 6.2 批处理方式

把 WT 和每个候选分别写为 `inputs/variants.fasta` 的独立记录。记录名应明确背景和替换，例如 `WT_<mutation>` 或 `N2V8_<mutation>`；替换标记中的原氨基酸必须与工作序列一致。

安装完成后，在 WSL2/Linux shell 中运行：

```bash
colabfold_batch --help
colabfold_batch inputs/variants.fasta results/colabfold \
  --model-type alphafold2_ptm \
  --num-models 5 \
  --num-recycle 3 \
  --random-seed 0 \
  --num-seeds 3
```

参数已按官方 CLI 核对；数值是本项目起始设置，并非活性优化参数。完整记录工具版本、权重、MSA 来源、模板设置和日志。[CLI 源码](https://github.com/sokrypton/ColabFold/blob/main/colabfold/batch.py)

单点突变比较尽量保持同一同源序列集合和模板策略；如复用 WT MSA，需要正确修改突变体 query 行，而不是把 WT query 原样提交。MSA 不同也可能造成模型差异。

## 7. 用 AlphaFold 3 分析 RNA 复合物

对 WT 和候选分别建立蛋白单体、蛋白–crRNA 二元复合物，以及蛋白–crRNA–匹配靶 RNA 三元复合物。使用实际实验 crRNA 的完整 repeat 和 spacer；靶 RNA 使用与实验一致的片段，记录其配对关系和边界。不同截短方式会影响比较。

网页方式：进入 [AlphaFold Server](https://alphafoldserver.com/)，分别添加 protein、RNA、RNA 实体，输入各自真实序列，保存完整结果。三元复合物更适合检查 RNA 识别与构象耦联，但静态模型不能证明非靶 RNA 的催化进入路径。[官方项目说明](https://github.com/google-deepmind/alphafold3)

本地方式：按官方安装文档部署 GPU 环境，并参考以下输入骨架。**三个占位符必须替换为真实序列后才能运行。** 该 JSON 是本地 `alphafold3` dialect，不能直接假定兼容网页服务器导入格式。

```json
{
  "name": "RfxCas13d_WT_ternary",
  "modelSeeds": [1, 2, 3],
  "sequences": [
    {"protein": {"id": "A", "sequence": "REPLACE_WITH_PROTEIN_SEQUENCE"}},
    {"rna": {"id": "B", "sequence": "REPLACE_WITH_CRRNA_SEQUENCE"}},
    {"rna": {"id": "C", "sequence": "REPLACE_WITH_TARGET_RNA_SEQUENCE"}}
  ],
  "dialect": "alphafold3",
  "version": 1
}
```

环境部署后，在 AlphaFold 3 仓库目录运行以下命令；路径需改为实际位置：

```bash
python run_alphafold.py \
  --json_path=/absolute/path/inputs/WT_ternary.json \
  --model_dir=/absolute/path/model_parameters \
  --db_dir=/absolute/path/databases \
  --output_dir=/absolute/path/results/af3
```

输入格式及运行方式依据：[输入文档](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)、[官方运行说明](https://github.com/google-deepmind/alphafold3)。二元模型删除靶 RNA 实体；突变体更换蛋白序列，并保持 RNA 与比较设置一致。

## 8. 分析模型与突变效应

1. **先检查可信度。** 查看突变附近 pLDDT、域间 PAE、RNA 界面置信度及跨模型一致性。低置信度区域只记为不确定，不据此断言发生结构重排。[AF3 输出说明](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)
2. **按域比对。** 在 PyMOL/ChimeraX 中先对齐稳定核心，再分别对齐 HEPN 等相关域。全蛋白 RMSD 可能掩盖域间运动，也可能把域运动误判为局部变形。
3. **检查局部变化。** 比较侧链冲突、氢键/盐桥、溶剂可及性和 RNA 接触。对比多个 WT 模型之间的差异与 WT–突变体差异；随机模型差异可能大于突变效应。
4. **补充稳定性评估。** 使用 FoldX 在相同参考骨架上建模，检查配对 WT 与突变体的能量差。它主要评估局部适配，不能替代激活动力学。

FoldX 最小示例（在准备好文件的目录运行）：

```bash
foldx --command=RepairPDB --pdb=reference.pdb
foldx --command=BuildModel --pdb=reference_Repair.pdb \
  --mutant-file=individual_list.txt --numberOfRuns=5
```

`individual_list.txt` 每行采用“原氨基酸 + PDB 链 + PDB 残基号 + 新氨基酸 + 分号”。例如 `YA142A;` 只表示**链 A 的 142 位确实为 Tyr 时**构建 Y→A；这是格式示例，不是本项目增强候选。组合替换用逗号分隔；不得未经映射直接使用 WT 编号。[突变文件格式](https://foldxsuite.crg.eu/parameter/mutant-file)

读取 `Dif_*` 中 mutant 相对配对 WT 的 ΔΔG，单位 kcal/mol；一般正值提示去稳定化。不要按单个能量数字判定增强；FoldX 对缺失坐标、低质量侧链和结构状态敏感。[BuildModel 文档](https://foldxsuite.crg.eu/command/BuildModel)

上述 FoldX 示例用于蛋白结构的稳定性初筛。若处理 RNA 复合物，需另行核对当前版本的 RNA 支持及参数，不能将删除 RNA 后的蛋白能量解释成 RNA 结合能。

## 9. 排序与实验验证

不设“pLDDT 越高越增强”或“ΔΔG 越负越增强”的排序规则。优先级取决于明确的机制依据、可靠的残基映射、跨模型一致性、局部结构是否合理及稳定性风险。结构无明显差异也不能排除功能变化。

验证时至少比较 WT、候选、催化失活对照、无靶 RNA 和不匹配/非靶向 crRNA 条件。分别测量：

- 匹配靶 RNA 的切割或敲低能力。
- 靶 RNA 激活后，独立非靶 RNA 报告底物的时间曲线和初始速率。
- 无靶背景切割、蛋白表达/可溶性、RNA 结合及细胞活性（适用于细胞实验）。

初始速率应在相同条件和线性区间比较，尽可能归一化到有效酶浓度。动力学参数只有在底物系列和相应模型适用时才能拟合；荧光 RFU/min 不能直接称为 kcat。细胞报告结果可作为表型证据，再用纯化蛋白实验区分表达效应与本征催化效应。

建议记录表：

| 候选/背景 | 位点映射 | 机制依据 | 局部置信度 | FoldX ΔΔG | 靶向活性 | 激活后旁切速率/WT | 无靶背景 | 结论 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| WT | 已确认 | 基线 | 待计算 | 参照 | 待测 | 1（测量后归一化） | 待测 | 基线 |
| 候选 | 待确认 | 待填写 | 待计算 | 待计算 | 待测 | 待测 | 待测 | 未验证 |

## 10. 建议输出与最短实施路线

以下目录是后续工作建议，本次仅创建本 README：

```text
inputs/       # WT/背景序列、候选序列、crRNA、靶 RNA、编号映射
structures/   # 原始 mmCIF、验证报告、参考状态说明
results/      # ColabFold、AF3、FoldX 的完整输出
analysis/     # 候选表、结构比较、截图与汇总
experiments/  # 原始时间曲线、归一化方式和重复实验结果
```

实际执行顺序：**确认 WT 与高保真背景 → 核对结构链和编号 → 提出少量单点假设 → 预测 WT 与候选 → 检查 RNA 复合物和稳定性 → 在统一条件下验证旁切与靶向活性 → 根据实测数据设计下一轮。**

可以报告“结构分析支持该候选进入实验验证”；只有获得对照充分的功能数据后，才能报告“该突变增强了 RfxCas13d collateral activity”。
