import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from datetime import datetime

# --------------------------
# 1. 配置参数（适配工作流输入）
# --------------------------
# 从命令行获取Excel文件路径（默认根目录的PDB_Dataset_Info_Full.xlsx）
if len(sys.argv) > 1:
    EXCEL_PATH = sys.argv[1]
else:
    EXCEL_PATH = "PDB_Dataset_Info_Full.xlsx"

# 结果输出目录（与工作流上传路径一致）
OUTPUT_DIR = "target_filter_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)  # 确保目录存在


# --------------------------
# 2. 核心筛选规则（抗细菌/抗病毒/抗虫）
# --------------------------
def define_filter_keywords():
    """定义三类靶点的筛选关键词（大小写不敏感）"""
    return {
        "抗细菌相关": [
            "bacterial", "bacteria", "e. coli", "staphylococcus", "streptococcus",
            "mycobacterium", "ribosome", "ribosomal rna", "bacillus", "salmonella",
            "enterococcus", "klebsiella", "pseudomonas"
        ],
        "抗病毒相关": [
            "viral", "virus", "hiv", "influenza", "hcv", "hbv", "sars",
            "coronavirus", "viral rna", "viral ribozyme", "zika", "dengue",
            "cmv", "ebv", "polio", "west nile"
        ],
        "抗虫相关": [
            "parasite", "malaria", "helminth", "trypanosoma", "leishmania",
            "schistosoma", "parasitic rna", "plasmodium", "toxoplasma", "cryptosporidium"
        ]
    }


def classify_target_type(desc_text, keywords_dict):
    """根据Description列文本，分类靶点类型"""
    desc_lower = str(desc_text).lower().strip()
    matched_types = []
    for target_type, keywords in keywords_dict.items():
        if any(keyword in desc_lower for keyword in keywords):
            matched_types.append(target_type)
    
    if not matched_types:
        return "非相关靶点"
    elif len(matched_types) > 1:
        return "多类别关联（" + " + ".join(matched_types) + "）"
    else:
        return matched_types[0]


# --------------------------
# 3. 主筛选流程
# --------------------------
def main():
    print(f"=== 开始RNA抗微生物靶点筛选 ===")
    print(f"Excel数据文件路径：{EXCEL_PATH}")
    print(f"筛选时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # 步骤1：读取PDB数据
    try:
        df = pd.read_excel(EXCEL_PATH)
        print(f"✅ 成功读取数据：共 {len(df)} 个PDB结构")
    except FileNotFoundError:
        print(f"❌ 错误：未找到文件 {EXCEL_PATH}，请检查路径是否正确")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 读取文件失败：{str(e)}")
        sys.exit(1)

    # 步骤2：筛选靶点
    keywords = define_filter_keywords()
    df["靶点关联类型"] = df["Description (描述)"].apply(lambda x: classify_target_type(x, keywords))
    
    # 保留仅抗细菌/抗病毒/抗虫相关的靶点（去重PDB ID）
    target_types_to_keep = ["抗细菌相关", "抗病毒相关", "抗虫相关", "多类别关联（抗细菌相关 + 抗病毒相关）", "多类别关联（抗细菌相关 + 抗虫相关）", "多类别关联（抗病毒相关 + 抗虫相关）"]
    filtered_df = df[df["靶点关联类型"].isin(target_types_to_keep)].drop_duplicates(subset=["PDB ID"], keep="first")
    print(f"✅ 筛选完成：共 {len(filtered_df)} 个抗微生物相关RNA靶点")

    # 步骤3：统计各类靶点数量
    target_count = filtered_df["靶点关联类型"].value_counts()
    print(f"\n📊 靶点分类统计：")
    for target_type, count in target_count.items():
        print(f"   - {target_type}：{count} 个")

    # 步骤4：生成结果文件
    # 4.1 总表（Excel，含所有信息）✅ 移除encoding参数，使用openpyxl引擎
    total_output_path = os.path.join(OUTPUT_DIR, "PDB_抗微生物RNA靶点总表.xlsx")
    filtered_df.to_excel(total_output_path, index=False, engine="openpyxl")
    
    # 4.2 分类子表（CSV，便于单独分析）
    for target_type in ["抗细菌相关", "抗病毒相关", "抗虫相关"]:
        sub_df = filtered_df[filtered_df["靶点关联类型"] == target_type]
        sub_output_path = os.path.join(OUTPUT_DIR, f"PDB_{target_type}RNA靶点表.csv")
        sub_df.to_csv(sub_output_path, index=False, encoding="utf-8-sig")
    
    # 4.3 可视化图表（靶点数量分布）
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']  # 避免中文乱码
    plt.figure(figsize=(12, 6))
    sns.countplot(
        x="靶点关联类型", 
        data=filtered_df, 
        palette=["#2ecc71", "#3498db", "#e74c3c", "#f39c12", "#9b59b6", "#1abc9c"]
    )
    plt.title("Distribution of Antimicrobial-Related RNA Targets in PDB", fontsize=14, pad=20)
    plt.xlabel("Target Type", fontsize=12)
    plt.ylabel("Number of PDB Structures", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    # 添加数值标签
    for i, v in enumerate(filtered_df["靶点关联类型"].value_counts().values):
        plt.text(i, v + 0.3, str(v), ha="center", fontsize=11, fontweight="bold")
    plt.tight_layout()
    plot_output_path = os.path.join(OUTPUT_DIR, "RNA靶点分类统计图.png")
    plt.savefig(plot_output_path, dpi=300, bbox_inches="tight")
    plt.close()

    # 4.4 筛选报告（TXT，含详细说明）
    report_content = f"""
# RNA抗微生物靶点筛选报告
=====================================
筛选时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
数据来源：{EXCEL_PATH}
原始数据量：{len(df)} 个PDB结构
筛选后数据量：{len(filtered_df)} 个抗微生物相关RNA靶点

## 一、筛选规则
1. 抗细菌相关靶点：匹配关键词（bacterial/ribosome/E. coli等）
2. 抗病毒相关靶点：匹配关键词（viral/HIV/influenza等）
3. 抗虫相关靶点：匹配关键词（parasite/malaria/plasmodium等）
4. 去重规则：按PDB ID去重，保留每个靶点的唯一记录

## 二、靶点分类统计
{target_count.to_string()}

## 三、核心靶点示例
### 1. 抗细菌相关
{filtered_df[filtered_df["靶点关联类型"] == "抗细菌相关"][["PDB ID", "Description (描述)", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

### 2. 抗病毒相关
{filtered_df[filtered_df["靶点关联类型"] == "抗病毒相关"][["PDB ID", "Description (描述)", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

### 3. 抗虫相关
{filtered_df[filtered_df["靶点关联类型"] == "抗虫相关"][["PDB ID", "Description (描述)", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

## 四、使用建议
1. 优先研究“抗细菌相关”靶点（数量最多，成药性高，如细菌核糖体rRNA）
2. 抗病毒靶点中，HIV-1 TAR RNA、流感病毒RNA可重点关注（临床需求大）
3. 抗虫靶点数量较少，可补充PDB数据库搜索（关键词：parasite RNA / plasmodium RNA）
4. 结果文件说明：
   - 总表：含所有筛选靶点的完整信息（Excel）
   - 分类表：按类型拆分的CSV文件，便于单独分析
   - 统计图：靶点分布可视化，用于汇报展示
"""
    report_output_path = os.path.join(OUTPUT_DIR, "RNA靶点筛选报告.txt")
    with open(report_output_path, "w", encoding="utf-8") as f:
        f.write(report_content)

    print(f"\n📄 所有结果已保存到 {OUTPUT_DIR} 目录：")
    print(f"   - 总表：PDB_抗微生物RNA靶点总表.xlsx")
    print(f"   - 分类表：PDB_抗细菌/抗病毒/抗虫相关RNA靶点表.csv")
    print(f"   - 图表：RNA靶点分类统计图.png")
    print(f"   - 报告：RNA靶点筛选报告.txt")
    print(f"\n=== 筛选流程全部完成 ===")


if __name__ == "__main__":
    main()
