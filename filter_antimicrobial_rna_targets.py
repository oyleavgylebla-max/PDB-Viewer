import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from datetime import datetime
from deep_translator import GoogleTranslator  # ✅ 新增：翻译库
from time import sleep  # ✅ 新增：避免翻译请求过快

# --------------------------
# 1. 配置参数（适配工作流输入）
# --------------------------
if len(sys.argv) > 1:
    EXCEL_PATH = sys.argv[1]
else:
    EXCEL_PATH = "PDB_Dataset_Info_Full.xlsx"

OUTPUT_DIR = "target_filter_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 初始化翻译器（设置超时和重试）
translator = GoogleTranslator(source='en', target='zh-CN')


# --------------------------
# 2. 新增功能1：提取RNA靶点名称
# --------------------------
def extract_rna_target_name(desc_text):
    """从英文Description中智能提取RNA靶点名称"""
    desc_lower = str(desc_text).lower().strip()
    
    # 优先级1：核糖开关（Riboswitch）
    if "riboswitch" in desc_lower:
        # 提取具体类型（如THF riboswitch → THF核糖开关）
        riboswitch_match = re.search(r'(\w+)\s+riboswitch', desc_lower)
        if riboswitch_match:
            return f"{riboswitch_match.group(1).upper()}核糖开关"
        else:
            return "细菌核糖开关"
    
    # 优先级2：核糖体RNA（Ribosomal RNA / rRNA）
    elif any(kw in desc_lower for kw in ["ribosomal", "ribosome", "rrna"]):
        # 提取物种信息（如E. coli → 大肠杆菌核糖体RNA）
        if "e. coli" in desc_lower:
            return "大肠杆菌核糖体RNA"
        elif "bacillus" in desc_lower:
            return "芽孢杆菌核糖体RNA"
        elif "plasmodium" in desc_lower:
            return "疟原虫核糖体RNA"
        elif "leishmania" in desc_lower:
            return "利什曼原虫核糖体RNA"
        else:
            return "核糖体RNA"
    
    # 优先级3：RNA适配体（Aptamer）
    elif "aptamer" in desc_lower:
        # 提取配体信息（如GTP aptamer → GTP RNA适配体）
        aptamer_match = re.search(r'(\w+)\s+aptamer', desc_lower)
        if aptamer_match:
            return f"{aptamer_match.group(1).upper()} RNA适配体"
        else:
            return "RNA适配体"
    
    # 优先级4：病毒RNA（Viral RNA）
    elif any(kw in desc_lower for kw in ["hiv", "influenza", "hcv", "viral rna"]):
        if "hiv" in desc_lower:
            return "HIV-1 TAR RNA"
        elif "influenza" in desc_lower:
            return "流感病毒RNA"
        elif "hcv" in desc_lower:
            return "丙肝病毒IRES RNA"
        else:
            return "病毒RNA"
    
    # 优先级5：其他特殊结构
    elif "g-quadruplex" in desc_lower or "g4" in desc_lower:
        return "G-四联体RNA"
    elif "ribozyme" in desc_lower:
        return "RNA核酶"
    else:
        return "未知RNA靶点"


# --------------------------
# 3. 新增功能2：英文描述翻译成中文（带重试机制）
# --------------------------
def translate_description(eng_text, max_retries=3):
    """将英文Description翻译成中文，避免请求过快"""
    if pd.isna(eng_text) or str(eng_text).strip() == "":
        return "无描述"
    
    for attempt in range(max_retries):
        try:
            # 限制翻译长度（Google Translate单次最多5000字符）
            text_to_translate = str(eng_text)[:4500]
            zh_text = translator.translate(text_to_translate)
            sleep(0.5)  # 每次翻译间隔0.5秒，避免被限流
            return zh_text
        except Exception as e:
            print(f"   ⚠️  翻译失败（尝试 {attempt+1}/{max_retries}）：{str(e)[:50]}")
            sleep(1)  # 重试前等待1秒
    
    # 多次失败后返回原文
    return str(eng_text)


# --------------------------
# 4. 核心筛选规则（抗细菌/抗病毒/抗虫）
# --------------------------
def define_filter_keywords():
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
# 5. 主筛选流程（升级完整版）
# --------------------------
def main():
    print(f"=== 开始RNA抗微生物靶点筛选（升级版：含靶点名称+中文描述） ===")
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
    
    target_types_to_keep = ["抗细菌相关", "抗病毒相关", "抗虫相关", "多类别关联（抗细菌相关 + 抗病毒相关）", "多类别关联（抗细菌相关 + 抗虫相关）", "多类别关联（抗病毒相关 + 抗虫相关）"]
    filtered_df = df[df["靶点关联类型"].isin(target_types_to_keep)].drop_duplicates(subset=["PDB ID"], keep="first").copy()
    print(f"✅ 筛选完成：共 {len(filtered_df)} 个抗微生物相关RNA靶点")

    # 步骤3：✨ 新增：提取RNA靶点名称
    print(f"\n🔍 正在提取RNA靶点名称...")
    filtered_df["RNA靶点名称"] = filtered_df["Description (描述)"].apply(extract_rna_target_name)
    print(f"   ✅ 靶点名称提取完成")

    # 步骤4：✨ 新增：翻译Description为中文
    print(f"\n🌍 正在翻译英文描述为中文（约需1-2分钟，请勿中断）...")
    filtered_df["描述（中文）"] = filtered_df["Description (描述)"].apply(translate_description)
    print(f"   ✅ 中文翻译完成")

    # 步骤5：统计各类靶点数量
    target_count = filtered_df["靶点关联类型"].value_counts()
    print(f"\n📊 靶点分类统计：")
    for target_type, count in target_count.items():
        print(f"   - {target_type}：{count} 个")

    # 步骤6：生成结果文件（调整列顺序，优先展示中文信息）
    # 重新排列列顺序：PDB ID → RNA靶点名称 → 靶点关联类型 → 描述（中文）→ 原始英文描述 → 其他列
    output_cols = ["PDB ID", "RNA靶点名称", "靶点关联类型", "描述（中文）", "Description (描述)", "Ligands (对应小分子)", "Publication (文章出处)"]
    # 确保所有列都存在（避免因原始数据列名不同报错）
    output_cols = [col for col in output_cols if col in filtered_df.columns]
    filtered_df = filtered_df[output_cols]

    # 6.1 总表（Excel）
    total_output_path = os.path.join(OUTPUT_DIR, "PDB_抗微生物RNA靶点总表_中文版.xlsx")
    filtered_df.to_excel(total_output_path, index=False, engine="openpyxl")
    
    # 6.2 分类子表（CSV）
    for target_type in ["抗细菌相关", "抗病毒相关", "抗虫相关"]:
        sub_df = filtered_df[filtered_df["靶点关联类型"] == target_type]
        sub_output_path = os.path.join(OUTPUT_DIR, f"PDB_{target_type}RNA靶点表_中文版.csv")
        sub_df.to_csv(sub_output_path, index=False, encoding="utf-8-sig")
    
    # 6.3 可视化图表（按RNA靶点名称统计）
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.figure(figsize=(14, 7))
    # 统计Top 10 RNA靶点名称
    top_targets = filtered_df["RNA靶点名称"].value_counts().head(10)
    sns.barplot(x=top_targets.values, y=top_targets.index, palette="viridis")
    plt.title("Top 10 Antimicrobial RNA Targets in PDB", fontsize=14, pad=20)
    plt.xlabel("Number of PDB Structures", fontsize=12)
    plt.ylabel("RNA Target Name", fontsize=12)
    # 添加数值标签
    for i, v in enumerate(top_targets.values):
        plt.text(v + 0.2, i, str(v), va="center", fontsize=11, fontweight="bold")
    plt.tight_layout()
    plot_output_path = os.path.join(OUTPUT_DIR, "RNA靶点名称分布统计图.png")
    plt.savefig(plot_output_path, dpi=300, bbox_inches="tight")
    plt.close()

    # 6.4 筛选报告（中文版）
    report_content = f"""
# RNA抗微生物靶点筛选报告（中文版）
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
5. 新增功能：自动提取RNA靶点名称 + 英文描述中文化

## 二、靶点分类统计
{target_count.to_string()}

## 三、Top 5 RNA靶点名称统计
{top_targets.head(5).to_string()}

## 四、核心靶点示例（含中文描述）
### 1. 抗细菌相关
{filtered_df[filtered_df["靶点关联类型"] == "抗细菌相关"][["PDB ID", "RNA靶点名称", "描述（中文）", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

### 2. 抗病毒相关
{filtered_df[filtered_df["靶点关联类型"] == "抗病毒相关"][["PDB ID", "RNA靶点名称", "描述（中文）", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

### 3. 抗虫相关
{filtered_df[filtered_df["靶点关联类型"] == "抗虫相关"][["PDB ID", "RNA靶点名称", "描述（中文）", "Ligands (对应小分子)"]].head(3).to_string(index=False)}

## 五、使用建议
1. 优先研究“核糖体RNA”和“核糖开关”（数量最多，成药性高）
2. 抗病毒靶点中，“HIV-1 TAR RNA”和“流感病毒RNA”可重点关注
3. 抗虫靶点数量较少，可补充PDB数据库搜索（关键词：parasite RNA）
4. 结果文件说明：
   - 总表：含所有筛选靶点的完整信息（Excel，中文版）
   - 分类表：按类型拆分的CSV文件，便于单独分析
   - 统计图：RNA靶点名称分布可视化，用于汇报展示
"""
    report_output_path = os.path.join(OUTPUT_DIR, "RNA靶点筛选报告_中文版.txt")
    with open(report_output_path, "w", encoding="utf-8") as f:
        f.write(report_content)

    print(f"\n📄 所有结果已保存到 {OUTPUT_DIR} 目录：")
    print(f"   - 总表：PDB_抗微生物RNA靶点总表_中文版.xlsx")
    print(f"   - 分类表：PDB_抗细菌/抗病毒/抗虫相关RNA靶点表_中文版.csv")
    print(f"   - 图表：RNA靶点名称分布统计图.png")
    print(f"   - 报告：RNA靶点筛选报告_中文版.txt")
    print(f"\n=== 筛选流程全部完成（升级版） ===")


if __name__ == "__main__":
    main()
