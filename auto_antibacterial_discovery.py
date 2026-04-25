import pandas as pd
import requests
import urllib.parse
import time
import json
import os
import sys
from functools import lru_cache
from datetime import datetime
# =============================================================================
# 🧬 全RNA靶点药物重定位自动化分析引擎
# 支持两种模式：全RNA分析 / 抗菌专项分析
# 完全兼容GitHub Actions v5工作流
# 2024优化版：全量适应症补全 + 官方靶点字段新增
# =============================================================================
# --------------------------
# 配置与缓存系统
# --------------------------
SMILES_CACHE_FILE = "smiles_cache.json"
DRUG_SEARCH_CACHE_FILE = "drug_search_cache.json"
# 输出目录（与GitHub Actions工作流完全一致）
ALL_RNA_OUTPUT_DIR = "all_rna_drug_discovery_results"
ANTIBACTERIAL_OUTPUT_DIR = "antibacterial_results"
# ATC分类与适应症映射（完整版）
ATC_FULL_MAP = {
    "A01": "口腔病用药", "A02": "治疗胃酸相关疾病的药物", "A03": "治疗功能性胃肠道疾病的药物",
    "A04": "止吐药和止恶心药", "A05": "胆和肝治疗药", "A06": "轻泻药", "A07": "肠道抗感染药和肠道消炎药",
    "A08": "减肥药", "A09": "消化药，包括酶", "A10": "糖尿病用药", "A11": "维生素类", "A12": "矿物质补充剂",
    "A13": "滋补药", "A14": "全身用蛋白同化类固醇", "A15": "食欲刺激药", "A16": "其他消化道和代谢药物",
    "B01": "抗血栓药", "B02": "止血药", "B03": "抗贫血药", "B05": "血液代用品和灌注液",
    "B06": "其他血液系统用药",
    "C01": "心脏治疗药", "C02": "抗高血压药", "C03": "利尿药", "C04": "外周血管扩张药",
    "C05": "血管保护剂", "C07": "β受体阻滞剂", "C08": "钙通道阻滞剂", "C09": "作用于肾素-血管紧张素系统的药物",
    "C10": "血脂调节剂",
    "D01": "皮肤用抗真菌药", "D02": "润肤剂和保护剂", "D03": "皮肤用皮质类固醇",
    "D04": "止痒药，包括抗组胺药、麻醉药", "D05": "银屑病用药", "D06": "皮肤用抗生素和化疗药",
    "D07": "皮肤用皮质类固醇和抗生素的复方制剂", "D08": "皮肤用消毒剂和防腐剂",
    "D09": "伤口敷料和保护剂", "D10": "痤疮用药", "D11": "其他皮肤科用药",
    "J01": "全身用抗菌药", "J02": "全身用抗真菌药", "J04": "抗分枝杆菌药", "J05": "全身用抗病毒药",
    "J06": "免疫血清和免疫球蛋白", "J07": "疫苗",
    "L01": "抗肿瘤药", "L02": "内分泌治疗药", "L03": "免疫刺激剂", "L04": "免疫抑制剂",
    "M01": "抗炎和抗风湿药", "M02": "局部用肌肉骨骼系统药物", "M03": "肌肉松弛药",
    "M04": "抗痛风药", "M05": "治疗骨病的药物",
    "N01": "麻醉药", "N02": "镇痛药", "N03": "抗癫痫药", "N04": "抗帕金森病药",
    "N05": "精神安定药", "N06": "精神兴奋药", "N07": "其他神经系统药物",
    "R01": "鼻腔用药", "R02": "咽喉用药", "R03": "用于阻塞性气道疾病的药物",
    "R05": "咳嗽和感冒用药", "R06": "全身用抗组胺药", "R07": "其他呼吸系统药物"
}
ATC_CLASS_MAP = {
    "A": "消化系统及代谢药", "B": "血液和造血系统药物", "C": "心血管系统药物",
    "D": "皮肤科用药", "G": "泌尿生殖系统药和性激素", "H": "全身用激素类制剂",
    "J": "全身用抗感染药", "L": "抗肿瘤药和免疫调节剂", "M": "肌肉-骨骼系统药物",
    "N": "神经系统药物", "P": "抗寄生虫药", "R": "呼吸系统药物",
    "S": "感觉器官药物", "V": "其他药品"
}
INDICATION_CN_MAP = {
    "cancer": "癌症", "tumor": "肿瘤", "carcinoma": "恶性肿瘤", "leukemia": "白血病",
    "hiv": "艾滋病", "influenza": "流感", "hepatitis": "肝炎", "bacterial infection": "细菌感染",
    "diabetes": "糖尿病", "obesity": "肥胖症", "hypertension": "高血压", "cardiovascular": "心血管疾病",
    "alzheimer": "阿尔茨海默病", "parkinson": "帕金森病", "arthritis": "关节炎", "asthma": "哮喘",
    "depression": "抑郁症", "pain": "疼痛", "inflammation": "炎症", "fungal infection": "真菌感染",
    "virus infection": "病毒感染", "thromboembolism": "血栓栓塞", "stroke": "中风",
    "psychosis": "精神病", "schizophrenia": "精神分裂症", "epilepsy": "癫痫",
    "gastroesophageal reflux": "胃食管反流", "ulcer": "消化道溃疡", "anemia": "贫血",
    "glaucoma": "青光眼", "migraine": "偏头痛", "osteoporosis": "骨质疏松症",
    "infection": "感染性疾病", "autoimmune": "自身免疫性疾病", "neuropathic pain": "神经病理性疼痛",
    "tuberculosis": "结核病", "malaria": "疟疾", "pneumonia": "肺炎", "sepsis": "脓毒症"
}
# ====================== 【核心优化1：全量药物适应症字典】======================
# 覆盖所有匹配药物+盐型/衍生物，彻底解决「暂无明确适应症」问题
MANUAL_DRUG_INDICATION = {
    # 嘌呤/嘧啶类抗肿瘤药
    "THIOGUANINE": "急性髓系白血病、急性淋巴细胞白血病",
    "THIOGUANINE TABLETS": "急性髓系白血病、急性淋巴细胞白血病",
    "PEMETREXED": "恶性胸膜间皮瘤、非小细胞肺癌",
    "PEMETREXED DIPOTASSIUM": "恶性胸膜间皮瘤、非小细胞肺癌",
    "PEMETREXED MONOHYDRATE": "恶性胸膜间皮瘤、非小细胞肺癌",
    "PEMETREXED TROMETHAMINE": "恶性胸膜间皮瘤、非小细胞肺癌",
    "PEMETREXED DISODIUM": "恶性胸膜间皮瘤、非小细胞肺癌",
    "PEMETREXED DISODIUM HEMIPENTAHYDRATE": "恶性胸膜间皮瘤、非小细胞肺癌",
    "FLUDARABINE PHOSPHATE": "慢性淋巴细胞白血病、非霍奇金淋巴瘤",
    "CLADRIBINE": "毛细胞白血病、慢性淋巴细胞白血病、急性髓系白血病",
    "CLOFARABINE": "儿童急性淋巴细胞白血病",
    "PENTOSTATIN": "毛细胞白血病、慢性淋巴细胞白血病、皮肤T细胞淋巴瘤",
    "NELARABINE": "T细胞急性淋巴细胞白血病、T细胞淋巴母细胞淋巴瘤",
    "CYTARABINE": "急性髓系白血病、急性淋巴细胞白血病、非霍奇金淋巴瘤",
    "AZACITIDINE": "骨髓增生异常综合征、急性髓系白血病、慢性粒单核细胞白血病",
    "DECITABINE": "骨髓增生异常综合征、急性髓系白血病",
    "FLOXURIDINE": "肝癌、胃肠道恶性肿瘤的动脉灌注化疗",
    "TRIFLURIDINE": "单纯疱疹病毒性角膜炎、联合替匹嘧啶用于转移性结直肠癌",
    "MITOXANTRONE": "急性白血病、恶性淋巴瘤、乳腺癌、前列腺癌、多发性硬化",
    "MITOXANTRONE HYDROCHLORIDE": "急性白血病、恶性淋巴瘤、乳腺癌、前列腺癌、多发性硬化",
    # 氨基糖苷类抗菌药
    "PAROMOMYCIN": "肠道阿米巴病、细菌性痢疾、隐孢子虫病、肠道细菌感染",
    "PAROMOMYCIN SULFATE": "肠道阿米巴病、细菌性痢疾、隐孢子虫病、肠道细菌感染",
    "NETILMICIN": "革兰氏阴性菌引起的呼吸道感染、泌尿道感染、皮肤软组织感染、败血症",
    "NETILMICIN SULFATE": "革兰氏阴性菌引起的呼吸道感染、泌尿道感染、皮肤软组织感染、败血症",
    "KANAMYCIN": "敏感革兰氏阴性菌所致严重感染、结核病二线治疗",
    "KANAMYCIN SULFATE": "敏感革兰氏阴性菌所致严重感染、结核病二线治疗",
    "GENTAMICIN": "革兰氏阴性菌所致的严重感染",
    "GENTAMICIN SULFATE": "革兰氏阴性菌所致的严重感染",
    "TOBRAMYCIN": "铜绿假单胞菌等革兰氏阴性菌所致的呼吸道、泌尿道、皮肤软组织、腹腔感染及败血症",
    "TOBRAMYCIN SULFATE": "铜绿假单胞菌等革兰氏阴性菌所致的呼吸道、泌尿道、皮肤软组织、腹腔感染及败血症",
    "AMIKACIN": "耐药革兰氏阴性菌所致严重感染、败血症、肺炎、中枢神经系统感染",
    "AMIKACIN SULFATE": "耐药革兰氏阴性菌所致严重感染、败血症、肺炎、中枢神经系统感染",
    "STREPTOMYCIN": "结核病一线治疗、鼠疫、布鲁菌病、兔热病",
    "STREPTOMYCIN SULFATE": "结核病一线治疗、鼠疫、布鲁菌病、兔热病",
    "PLAZOMICIN": "多重耐药革兰氏阴性菌所致严重感染、败血症、医院获得性肺炎",
    "PLAZOMICIN SULFATE": "多重耐药革兰氏阴性菌所致严重感染、败血症、医院获得性肺炎",
    # 大环内酯类抗菌药
    "ERYTHROMYCIN": "呼吸道感染、皮肤软组织感染、支原体肺炎、衣原体感染、军团菌病",
    # 四环素类抗菌药
    "TETRACYCLINE": "立克次体病、支原体肺炎、衣原体感染、布鲁菌病、霍乱",
    "TETRACYCLINE HYDROCHLORIDE": "立克次体病、支原体肺炎、衣原体感染、布鲁菌病、霍乱",
    "TETRACYCLINE PHOSPHATE COMPLEX": "立克次体病、支原体肺炎、衣原体感染、布鲁菌病、霍乱",
    "SARECYCLINE": "社区获得性细菌性肺炎、急性细菌性皮肤和皮肤结构感染",
    "SARECYCLINE HYDROCHLORIDE": "社区获得性细菌性肺炎、急性细菌性皮肤和皮肤结构感染",
    "OMADACYCLINE": "社区获得性细菌性肺炎、急性细菌性皮肤和皮肤结构感染",
    "OMADACYCLINE TOSYLATE": "社区获得性细菌性肺炎、急性细菌性皮肤和皮肤结构感染",
    "ERAVACYCLINE": "复杂性腹腔内感染、复杂性尿路感染、急性细菌性皮肤感染",
    "ERAVACYCLINE DIHYDROCHLORIDE": "复杂性腹腔内感染、复杂性尿路感染、急性细菌性皮肤感染",
    "MINOCYCLINE": "中重度痤疮、非淋菌性尿道炎、皮肤软组织感染、脑膜炎球菌带菌者清除",
    "MINOCYCLINE HYDROCHLORIDE": "中重度痤疮、非淋菌性尿道炎、皮肤软组织感染、脑膜炎球菌带菌者清除",
    "DOXYCYCLINE": "痤疮、立克次体病、支原体/衣原体感染、霍乱、布鲁菌病、疟疾预防",
    "DOXYCYCLINE HYCLATE": "痤疮、立克次体病、支原体/衣原体感染、霍乱、布鲁菌病、疟疾预防",
    "DOXYCYCLINE CALCIUM": "痤疮、立克次体病、支原体/衣原体感染、霍乱、布鲁菌病、疟疾预防",
    "METHACYCLINE": "敏感菌所致呼吸道感染、尿路感染、皮肤软组织感染",
    "METHACYCLINE HYDROCHLORIDE": "敏感菌所致呼吸道感染、尿路感染、皮肤软组织感染",
    "OXYTETRACYCLINE": "敏感菌所致呼吸道感染、皮肤感染、立克次体病、衣原体感染",
    "OXYTETRACYCLINE HYDROCHLORIDE": "敏感菌所致呼吸道感染、皮肤感染、立克次体病、衣原体感染",
    "OXYTETRACYCLINE CALCIUM": "敏感菌所致呼吸道感染、皮肤感染、立克次体病、衣原体感染",
    "CHLORTETRACYCLINE": "皮肤浅表细菌感染、细菌性结膜炎、麦粒肿、细菌性阴道病",
    "CHLORTETRACYCLINE HYDROCHLORIDE": "皮肤浅表细菌感染、细菌性结膜炎、麦粒肿、细菌性阴道病",
    "DEMECLOCYCLINE": "抗利尿激素分泌异常综合征、敏感菌所致感染",
    "DEMECLOCYCLINE HYDROCHLORIDE": "抗利尿激素分泌异常综合征、敏感菌所致感染",
    "MECLOCYCLINE SULFOSALICYLATE": "敏感菌所致呼吸道、泌尿道、皮肤软组织感染",
    # 维生素类
    "GUANIDINE": "重症肌无力、神经源性膀胱",
    "GUANIDINE HYDROCHLORIDE": "重症肌无力、神经源性膀胱",
    "RIBOFLAVIN 5'-PHOSPHATE": "核黄素缺乏症、口角炎、舌炎、脂溢性皮炎",
    "RIBOFLAVIN 5'-PHOSPHATE SODIUM": "核黄素缺乏症、口角炎、舌炎、脂溢性皮炎",
    "RIBOFLAVIN": "维生素B2缺乏症、口角炎、唇炎、阴囊炎、脂溢性皮炎",
    "FOLIC ACID": "巨幼细胞性贫血、妊娠期预防胎儿神经管畸形、高同型半胱氨酸血症",
    "LEUCOVORIN": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药、巨幼细胞性贫血",
    "LEUCOVORIN CALCIUM": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药、巨幼细胞性贫血",
    "LEVOLEUCOVORIN": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药",
    "LEVOLEUCOVORIN CALCIUM": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药",
    "LEVOMEFOLIC ACID": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药",
    "LEVOMEFOLATE CALCIUM": "甲氨蝶呤中毒解救、结直肠癌化疗辅助用药",
    "THIAMINE": "维生素B1缺乏症、脚气病、神经炎、韦尼克脑病",
    "THIAMINE HYDROCHLORIDE": "维生素B1缺乏症、脚气病、神经炎、韦尼克脑病",
    "THIAMINE ION": "维生素B1缺乏症、脚气病、神经炎、韦尼克脑病",
    "NIACINAMIDE": "烟酸缺乏症、糙皮病、口炎、舌炎",
    "BIOTIN": "生物素缺乏症、脂溢性皮炎、脱发、先天性生物素酶缺乏症",
    # 抗病毒药
    "VIDARABINE": "单纯疱疹病毒性角膜炎、带状疱疹、水痘病毒感染",
    "IDOXURIDINE": "单纯疱疹病毒性角膜炎、牛痘病毒性眼病",
    "TELBIVUDINE": "慢性乙型肝炎病毒感染",
    "ENTECAVIR": "慢性乙型肝炎病毒感染",
    "ENTECAVIR ANHYDROUS": "慢性乙型肝炎病毒感染",
    # 心血管系统药物
    "ADENOSINE": "阵发性室上性心动过速、室上性心律失常诊断与治疗",
    "REGADENOSON": "心肌灌注显像药物、冠状动脉疾病诊断用药",
    "DIQUAFOSOL": "干眼症、角膜上皮损伤",
    "DIQUAFOSOL TETRASODIUM": "干眼症、角膜上皮损伤",
    "PRAZOSIN": "高血压、良性前列腺增生",
    "PRAZOSIN HYDROCHLORIDE": "高血压、良性前列腺增生",
    "TERAZOSIN": "良性前列腺增生、高血压",
    "TERAZOSIN HYDROCHLORIDE": "良性前列腺增生、高血压",
    "DOXAZOSIN": "高血压、良性前列腺增生",
    "DOXAZOSIN MESYLATE": "高血压、良性前列腺增生",
    # 喹诺酮类抗菌药
    "CIPROFLOXACIN": "呼吸道、泌尿生殖道、消化道、皮肤软组织、腹腔感染及败血症、伤寒",
    "CIPROFLOXACIN HYDROCHLORIDE": "呼吸道、泌尿生殖道、消化道、皮肤软组织、腹腔感染及败血症、伤寒",
    # 其他抗感染药
    "POLIHEXANIDE": "皮肤、黏膜、创面的消毒防腐，伤口感染预防与治疗",
    # 激素类/内分泌药物
    "TAMOXIFEN": "雌激素受体阳性的乳腺癌、卵巢癌",
    "TAMOXIFEN CITRATE": "雌激素受体阳性的乳腺癌、卵巢癌",
    # 神经系统药物
    "PROMAZINE": "精神分裂症、躁狂症、恶心呕吐、顽固性呃逆",
    "PROMAZINE HYDROCHLORIDE": "精神分裂症、躁狂症、恶心呕吐、顽固性呃逆",
    "CHLORPROMAZINE": "精神分裂症、躁狂症、呕吐、顽固性呃逆、人工冬眠",
    "CHLORPROMAZINE HYDROCHLORIDE": "精神分裂症、躁狂症、呕吐、顽固性呃逆、人工冬眠",
    "ACETOPHENAZINE": "精神分裂症、偏执型精神障碍",
    "ACETOPHENAZINE MALEATE": "精神分裂症、偏执型精神障碍",
    "TRIFLUPROMAZINE": "精神分裂症、恶心呕吐、顽固性呃逆",
    "TRIFLUPROMAZINE HYDROCHLORIDE": "精神分裂症、恶心呕吐、顽固性呃逆",
    "PIPERACETAZINE": "精神分裂症、神经症、焦虑症",
    "CARPHENAZINE": "精神分裂症、偏执型精神障碍",
    "CARPHENAZINE MALEATE": "精神分裂症、偏执型精神障碍",
    "PROPIOMAZINE": "精神分裂症、恶心呕吐、术前镇静",
    "PROPIOMAZINE HYDROCHLORIDE": "精神分裂症、恶心呕吐、术前镇静",
    "IMIPRAMINE": "抑郁症、小儿遗尿症、注意缺陷多动障碍",
    "IMIPRAMINE HYDROCHLORIDE": "抑郁症、小儿遗尿症、注意缺陷多动障碍",
    "IMIPRAMINE PAMOATE": "抑郁症、小儿遗尿症",
    "PROCHLORPERAZINE": "精神分裂症、恶心呕吐、眩晕症",
    "PROCHLORPERAZINE MALEATE": "精神分裂症、恶心呕吐、眩晕症",
    "PROCHLORPERAZINE EDISYLATE": "精神分裂症、恶心呕吐、眩晕症",
    "PERPHENAZINE": "精神分裂症、恶心呕吐、顽固性呃逆",
    "TRIFLUOPERAZINE": "精神分裂症、焦虑症、恶心呕吐",
    "TRIFLUOPERAZINE HYDROCHLORIDE": "精神分裂症、焦虑症、恶心呕吐",
    "TRIMEPRAZINE": "荨麻疹、湿疹、过敏性皮炎、瘙痒症、术前镇静",
    "TRIMEPRAZINE TARTRATE": "荨麻疹、湿疹、过敏性皮炎、瘙痒症、术前镇静",
    "METHYLPROMAZINE": "过敏性疾病、瘙痒症、术前镇静、恶心呕吐",
    "CLOMIPRAMINE": "抑郁症、强迫症、恐怖症、焦虑症",
    "CLOMIPRAMINE HYDROCHLORIDE": "抑郁症、强迫症、恐怖症、焦虑症",
    "THIETHYLPERAZINE": "恶心呕吐、精神分裂症、焦虑症",
    "THIETHYLPERAZINE MALEATE": "恶心呕吐、精神分裂症、焦虑症",
    "THIETHYLPERAZINE MALATE": "恶心呕吐、精神分裂症、焦虑症",
    "FLUPHENAZINE": "精神分裂症、偏执型精神障碍",
    "FLUPHENAZINE HYDROCHLORIDE": "精神分裂症、偏执型精神障碍",
    "LEVOMEPROMAZINE": "精神分裂症、躁狂症、顽固性疼痛、术前镇静",
    "LEVODOPA": "帕金森病、帕金森综合征",
    "FLUORODOPA": "帕金森病的诊断与辅助治疗",
    "FLUORODOPA F 18": "帕金森病的PET显像诊断",
    "AMINOCAPROIC ACID": "纤溶亢进所致的出血、手术出血的预防与治疗",
    # 呼吸系统药物
    "THEOPHYLLINE": "支气管哮喘、慢性阻塞性肺疾病、喘息型支气管炎",
    "THEOPHYLLINE ANHYDROUS": "支气管哮喘、慢性阻塞性肺疾病、喘息型支气管炎",
    "THEOPHYLLINE GLYCINATE": "支气管哮喘、慢性阻塞性肺疾病",
    "THEOPHYLLINE SODIUM GLYCINATE": "支气管哮喘、慢性阻塞性肺疾病",
    "AMINOPHYLLINE": "支气管哮喘、喘息型支气管炎、慢性阻塞性肺疾病、心源性肺水肿",
    "CAFFEINE": "中枢性呼吸抑制、新生儿呼吸暂停、提神醒脑",
    "CAFFEINE CITRATE": "新生儿呼吸暂停",
    # 消化系统/代谢药物
    "LACTULOSE": "慢性功能性便秘、肝性脑病",
    "ARGININE HYDROCHLORIDE": "肝性脑病、代谢性碱中毒、精氨酸缺乏症",
    "GLUTAMINE": "胃肠道黏膜损伤、放化疗所致的胃肠道反应、高分解代谢状态",
    "CYSTEINE": "支气管炎、咳痰困难、乙酰氨基酚中毒解救、皮肤黏膜保护",
    "CYSTEINE HYDROCHLORIDE": "支气管炎、咳痰困难、乙酰氨基酚中毒解救",
    "SELENOMETHIONINE": "硒缺乏症的补充、克山病、大骨节病辅助治疗",
    "SELENOMETHIONINE SE 75": "胰腺疾病的显像诊断",
    "TYROSINE": "酪氨酸缺乏症、抑郁症、注意力缺陷障碍辅助治疗",
    "ETELCALCETIDE": "慢性肾病透析患者的继发性甲状旁腺功能亢进症",
    "ETELCALCETIDE HYDROCHLORIDE": "慢性肾病透析患者的继发性甲状旁腺功能亢进症",
    "ACETOHYDROXAMIC ACID": "尿路感染、尿路结石、尿素分解菌所致的感染",
    "AZELAIC ACID": "寻常痤疮、玫瑰痤疮、色素沉着过度",
    "POTASSIUM OXYBATE": "发作性睡病、猝倒症",
    "SODIUM OXYBATE": "发作性睡病、猝倒症",
    "OXYBATE": "发作性睡病、猝倒症",
    "MAGNESIUM OXYBATE": "发作性睡病辅助治疗",
    "CALCIUM OXYBATE": "发作性睡病辅助治疗",
    "AMINOLEVULINIC ACID": "光动力治疗，用于痤疮、皮肤癌、膀胱癌、尖锐湿疣",
    "AMINOLEVULINIC ACID HYDROCHLORIDE": "光动力治疗，用于痤疮、皮肤癌、膀胱癌",
    # 血液系统药物
    "TRIENTINE": "肝豆状核变性（威尔逊病），用于不能耐受青霉胺的患者",
    "TRIENTINE TETRAHYDROCHLORIDE": "肝豆状核变性（威尔逊病）",
    "TRIENTINE HYDROCHLORIDE": "肝豆状核变性（威尔逊病）",
    "AMIFOSTINE": "肿瘤放化疗的正常细胞保护剂，减少肾毒性、神经毒性",
    "DEXTRIP": "右旋糖酐，血容量补充、休克、血栓性疾病",
    "DEXTROSE": "低血糖症、能量补充、脱水、高钾血症",
    # 电解质/补充剂
    "SODIUM ACETATE": "代谢性酸中毒、电解质补充、血液透析缓冲液",
    "POTASSIUM ACETATE": "低钾血症的预防与治疗、电解质补充",
    "CALCIUM ACETATE": "慢性肾病所致的高磷血症、钙补充",
    "MAGNESIUM ACETATE": "镁缺乏症的补充、电解质紊乱",
    "ZINC ACETATE": "锌缺乏症、威尔逊病、腹泻辅助治疗",
    "ALUMINUM ACETATE": "皮肤湿疹、皮炎、溃疡的湿敷，收敛消炎",
    # 消毒/防腐药
    "BENZALKONIUM CHLORIDE": "皮肤、黏膜、器械的消毒防腐，创面感染预防",
    "CHLORHEXIDINE": "皮肤消毒、口腔感染、创面消毒、器械消毒",
    # 其他
    "KRYPTON": "肺部通气/灌注显像诊断、视网膜激光光凝治疗",
    "KRYPTON KR 81M": "肺部通气显像诊断",
    "SAPROPTERIN": "苯丙酮尿症、四氢生物蝶呤缺乏症",
    "SAPROPTERIN DIHYDROCHLORIDE": "苯丙酮尿症、四氢生物蝶呤缺乏症",
    "LISDEXAMFETAMINE": "注意缺陷多动障碍、暴食症",
    "LISDEXAMFETAMINE DIMESYLATE": "注意缺陷多动障碍、暴食症",
    "TEMOPORFIN": "光动力治疗，用于晚期头颈部癌、食管癌、肺癌",
}

# ====================== 【核心优化2：全量药物官方靶点字典】======================
# 与适应症字典一一对应，覆盖所有药物+盐型，官方药典/说明书靶点
MANUAL_DRUG_TARGET = {
    # 嘌呤/嘧啶类抗肿瘤药
    "THIOGUANINE": "嘌呤核苷磷酸化酶、DNA合成酶、嘌呤代谢相关酶",
    "THIOGUANINE TABLETS": "嘌呤核苷磷酸化酶、DNA合成酶、嘌呤代谢相关酶",
    "PEMETREXED": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "PEMETREXED DIPOTASSIUM": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "PEMETREXED MONOHYDRATE": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "PEMETREXED TROMETHAMINE": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "PEMETREXED DISODIUM": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "PEMETREXED DISODIUM HEMIPENTAHYDRATE": "胸苷酸合成酶(TS)、二氢叶酸还原酶(DHFR)、甘氨酰胺核糖核苷酸甲酰转移酶(GARFT)",
    "FLUDARABINE PHOSPHATE": "DNA聚合酶、核糖核苷酸还原酶",
    "CLADRIBINE": "腺苷脱氨酶、DNA合成酶、核糖核苷酸还原酶",
    "CLOFARABINE": "核糖核苷酸还原酶、DNA聚合酶",
    "PENTOSTATIN": "腺苷脱氨酶(ADA)",
    "NELARABINE": "DNA合成酶、嘌呤核苷磷酸化酶",
    "CYTARABINE": "DNA聚合酶、核糖核苷酸还原酶",
    "AZACITIDINE": "DNA甲基转移酶(DNMT)",
    "DECITABINE": "DNA甲基转移酶(DNMT)",
    "FLOXURIDINE": "胸苷酸合成酶(TS)",
    "TRIFLURIDINE": "病毒DNA聚合酶、胸苷酸合成酶(TS)",
    "MITOXANTRONE": "拓扑异构酶II、DNA嵌入剂",
    "MITOXANTRONE HYDROCHLORIDE": "拓扑异构酶II、DNA嵌入剂",
    # 氨基糖苷类抗菌药
    "PAROMOMYCIN": "细菌30S核糖体亚基，抑制细菌蛋白质合成",
    "PAROMOMYCIN SULFATE": "细菌30S核糖体亚基，抑制细菌蛋白质合成",
    "NETILMICIN": "细菌30S核糖体亚基",
    "NETILMICIN SULFATE": "细菌30S核糖体亚基",
    "KANAMYCIN": "细菌16S rRNA、30S核糖体亚基",
    "KANAMYCIN SULFATE": "细菌16S rRNA、30S核糖体亚基",
    "GENTAMICIN": "细菌30S核糖体亚基",
    "GENTAMICIN SULFATE": "细菌30S核糖体亚基",
    "TOBRAMYCIN": "细菌30S核糖体亚基、16S rRNA",
    "TOBRAMYCIN SULFATE": "细菌30S核糖体亚基、16S rRNA",
    "AMIKACIN": "细菌30S核糖体亚基、16S rRNA",
    "AMIKACIN SULFATE": "细菌30S核糖体亚基、16S rRNA",
    "STREPTOMYCIN": "细菌16S rRNA、30S核糖体亚基",
    "STREPTOMYCIN SULFATE": "细菌16S rRNA、30S核糖体亚基",
    "PLAZOMICIN": "细菌16S rRNA、30S核糖体亚基",
    "PLAZOMICIN SULFATE": "细菌16S rRNA、30S核糖体亚基",
    # 大环内酯类抗菌药
    "ERYTHROMYCIN": "细菌50S核糖体亚基，抑制细菌蛋白质合成",
    # 四环素类抗菌药
    "TETRACYCLINE": "细菌30S核糖体亚基，抑制细菌蛋白质合成",
    "TETRACYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基，抑制细菌蛋白质合成",
    "TETRACYCLINE PHOSPHATE COMPLEX": "细菌30S核糖体亚基，抑制细菌蛋白质合成",
    "SARECYCLINE": "细菌30S核糖体亚基",
    "SARECYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "OMADACYCLINE": "细菌30S核糖体亚基",
    "OMADACYCLINE TOSYLATE": "细菌30S核糖体亚基",
    "ERAVACYCLINE": "细菌30S核糖体亚基",
    "ERAVACYCLINE DIHYDROCHLORIDE": "细菌30S核糖体亚基",
    "MINOCYCLINE": "细菌30S核糖体亚基",
    "MINOCYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "DOXYCYCLINE": "细菌30S核糖体亚基",
    "DOXYCYCLINE HYCLATE": "细菌30S核糖体亚基",
    "DOXYCYCLINE CALCIUM": "细菌30S核糖体亚基",
    "METHACYCLINE": "细菌30S核糖体亚基",
    "METHACYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "OXYTETRACYCLINE": "细菌30S核糖体亚基",
    "OXYTETRACYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "OXYTETRACYCLINE CALCIUM": "细菌30S核糖体亚基",
    "CHLORTETRACYCLINE": "细菌30S核糖体亚基",
    "CHLORTETRACYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "DEMECLOCYCLINE": "细菌30S核糖体亚基",
    "DEMECLOCYCLINE HYDROCHLORIDE": "细菌30S核糖体亚基",
    "MECLOCYCLINE SULFOSALICYLATE": "细菌30S核糖体亚基",
    # 维生素类
    "GUANIDINE": "烟碱型乙酰胆碱受体",
    "GUANIDINE HYDROCHLORIDE": "烟碱型乙酰胆碱受体",
    "RIBOFLAVIN 5'-PHOSPHATE": "黄素蛋白辅酶、黄素单核苷酸(FMN)合成靶点",
    "RIBOFLAVIN 5'-PHOSPHATE SODIUM": "黄素蛋白辅酶、黄素单核苷酸(FMN)合成靶点",
    "RIBOFLAVIN": "黄素单核苷酸(FMN)、黄素腺嘌呤二核苷酸(FAD)合成通路",
    "FOLIC ACID": "二氢叶酸还原酶、叶酸合成代谢通路",
    "LEUCOVORIN": "叶酸代谢通路、二氢叶酸还原酶",
    "LEUCOVORIN CALCIUM": "叶酸代谢通路、二氢叶酸还原酶",
    "LEVOLEUCOVORIN": "叶酸代谢通路、二氢叶酸还原酶",
    "LEVOLEUCOVORIN CALCIUM": "叶酸代谢通路、二氢叶酸还原酶",
    "LEVOMEFOLIC ACID": "叶酸代谢通路、二氢叶酸还原酶",
    "LEVOMEFOLATE CALCIUM": "叶酸代谢通路、二氢叶酸还原酶",
    "THIAMINE": "酮酸脱氢酶复合体、转酮醇酶，参与糖代谢",
    "THIAMINE HYDROCHLORIDE": "酮酸脱氢酶复合体、转酮醇酶，参与糖代谢",
    "THIAMINE ION": "酮酸脱氢酶复合体、转酮醇酶，参与糖代谢",
    "NIACINAMIDE": "烟酰胺腺嘌呤二核苷酸(NAD)、烟酰胺腺嘌呤二核苷酸磷酸(NADP)合成靶点",
    "BIOTIN": "乙酰辅酶A羧化酶、丙酮酸羧化酶等羧化酶的辅酶",
    # 抗病毒药
    "VIDARABINE": "病毒DNA聚合酶",
    "IDOXURIDINE": "病毒胸苷激酶、病毒DNA合成酶",
    "TELBIVUDINE": "乙型肝炎病毒逆转录酶",
    "ENTECAVIR": "乙型肝炎病毒逆转录酶",
    "ENTECAVIR ANHYDROUS": "乙型肝炎病毒逆转录酶",
    # 心血管系统药物
    "ADENOSINE": "腺苷A1受体、腺苷A2A受体",
    "REGADENOSON": "腺苷A2A受体",
    "DIQUAFOSOL": "P2Y2嘌呤能受体",
    "DIQUAFOSOL TETRASODIUM": "P2Y2嘌呤能受体",
    "PRAZOSIN": "α1肾上腺素能受体",
    "PRAZOSIN HYDROCHLORIDE": "α1肾上腺素能受体",
    "TERAZOSIN": "α1肾上腺素能受体",
    "TERAZOSIN HYDROCHLORIDE": "α1肾上腺素能受体",
    "DOXAZOSIN": "α1肾上腺素能受体",
    "DOXAZOSIN MESYLATE": "α1肾上腺素能受体",
    # 喹诺酮类抗菌药
    "CIPROFLOXACIN": "细菌DNA旋转酶（拓扑异构酶II）、拓扑异构酶IV",
    "CIPROFLOXACIN HYDROCHLORIDE": "细菌DNA旋转酶（拓扑异构酶II）、拓扑异构酶IV",
    # 其他抗感染药
    "POLIHEXANIDE": "细菌细胞膜、微生物核酸，破坏微生物结构完整性",
    # 激素类/内分泌药物
    "TAMOXIFEN": "雌激素受体(ER)，竞争性抑制雌激素结合",
    "TAMOXIFEN CITRATE": "雌激素受体(ER)，竞争性抑制雌激素结合",
    # 神经系统药物
    "PROMAZINE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "PROMAZINE HYDROCHLORIDE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "CHLORPROMAZINE": "多巴胺D2受体、5-HT受体、组胺H1受体、M胆碱受体、α肾上腺素能受体",
    "CHLORPROMAZINE HYDROCHLORIDE": "多巴胺D2受体、5-HT受体、组胺H1受体、M胆碱受体、α肾上腺素能受体",
    "ACETOPHENAZINE": "多巴胺D2受体",
    "ACETOPHENAZINE MALEATE": "多巴胺D2受体",
    "TRIFLUPROMAZINE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "TRIFLUPROMAZINE HYDROCHLORIDE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "PIPERACETAZINE": "多巴胺D2受体",
    "CARPHENAZINE": "多巴胺D2受体",
    "CARPHENAZINE MALEATE": "多巴胺D2受体",
    "PROPIOMAZINE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "PROPIOMAZINE HYDROCHLORIDE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "IMIPRAMINE": "5-羟色胺、去甲肾上腺素再摄取抑制剂",
    "IMIPRAMINE HYDROCHLORIDE": "5-羟色胺、去甲肾上腺素再摄取抑制剂",
    "IMIPRAMINE PAMOATE": "5-羟色胺、去甲肾上腺素再摄取抑制剂",
    "PROCHLORPERAZINE": "多巴胺D2受体、组胺H1受体、M胆碱受体",
    "PROCHLORPERAZINE MALEATE": "多巴胺D2受体、组胺H1受体、M胆碱受体",
    "PROCHLORPERAZINE EDISYLATE": "多巴胺D2受体、组胺H1受体、M胆碱受体",
    "PERPHENAZINE": "多巴胺D2受体、5-HT受体",
    "TRIFLUOPERAZINE": "多巴胺D2受体",
    "TRIFLUOPERAZINE HYDROCHLORIDE": "多巴胺D2受体",
    "TRIMEPRAZINE": "组胺H1受体、多巴胺D2受体、M胆碱受体",
    "TRIMEPRAZINE TARTRATE": "组胺H1受体、多巴胺D2受体、M胆碱受体",
    "METHYLPROMAZINE": "组胺H1受体、多巴胺D2受体",
    "CLOMIPRAMINE": "5-羟色胺再摄取强效抑制剂，去甲肾上腺素再摄取抑制剂",
    "CLOMIPRAMINE HYDROCHLORIDE": "5-羟色胺再摄取强效抑制剂，去甲肾上腺素再摄取抑制剂",
    "THIETHYLPERAZINE": "多巴胺D2受体、组胺H1受体",
    "THIETHYLPERAZINE MALEATE": "多巴胺D2受体、组胺H1受体",
    "THIETHYLPERAZINE MALATE": "多巴胺D2受体、组胺H1受体",
    "FLUPHENAZINE": "多巴胺D2受体",
    "FLUPHENAZINE HYDROCHLORIDE": "多巴胺D2受体",
    "LEVOMEPROMAZINE": "多巴胺D2受体、组胺H1受体、α肾上腺素能受体",
    "LEVODOPA": "多巴胺前体，经多巴脱羧酶转化为多巴胺",
    "FLUORODOPA": "多巴胺前体，芳香族L-氨基酸脱羧酶底物",
    "FLUORODOPA F 18": "芳香族L-氨基酸脱羧酶底物，PET显像剂",
    "AMINOCAPROIC ACID": "纤溶酶原激活剂抑制剂，抑制纤溶酶原激活",
    # 呼吸系统药物
    "THEOPHYLLINE": "磷酸二酯酶、腺苷受体，松弛支气管平滑肌",
    "THEOPHYLLINE ANHYDROUS": "磷酸二酯酶、腺苷受体，松弛支气管平滑肌",
    "THEOPHYLLINE GLYCINATE": "磷酸二酯酶、腺苷受体，松弛支气管平滑肌",
    "THEOPHYLLINE SODIUM GLYCINATE": "磷酸二酯酶、腺苷受体，松弛支气管平滑肌",
    "AMINOPHYLLINE": "磷酸二酯酶、腺苷受体，松弛支气管平滑肌",
    "CAFFEINE": "腺苷受体、磷酸二酯酶，中枢神经兴奋",
    "CAFFEINE CITRATE": "腺苷受体、磷酸二酯酶，中枢神经兴奋",
    # 消化系统/代谢药物
    "LACTULOSE": "肠道菌群代谢靶点，降低肠道pH，促进肠道蠕动",
    "ARGININE HYDROCHLORIDE": "精氨酸酶、一氧化氮合酶，促进尿素循环",
    "GLUTAMINE": "谷氨酰胺酶、肠道黏膜细胞代谢靶点，促进黏膜修复",
    "CYSTEINE": "黏液溶解靶点，还原黏蛋白二硫键；谷胱甘肽合成前体",
    "CYSTEINE HYDROCHLORIDE": "黏液溶解靶点，还原黏蛋白二硫键；谷胱甘肽合成前体",
    "SELENOMETHIONINE": "谷胱甘肽过氧化物酶、硫氧还蛋白还原酶的活性中心组成成分",
    "SELENOMETHIONINE SE 75": "谷胱甘肽过氧化物酶底物，PET显像剂",
    "TYROSINE": "酪氨酸羟化酶、多巴脱羧酶，多巴胺/肾上腺素合成前体",
    "ETELCALCETIDE": "钙敏感受体(CaSR)激动剂",
    "ETELCALCETIDE HYDROCHLORIDE": "钙敏感受体(CaSR)激动剂",
    "ACETOHYDROXAMIC ACID": "尿素酶抑制剂，抑制细菌尿素分解",
    "AZELAIC ACID": "5α-还原酶、酪氨酸酶抑制剂，抗菌抗炎，抑制毛囊角化",
    "POTASSIUM OXYBATE": "GABA-B受体、谷氨酸脱羧酶，中枢抑制",
    "SODIUM OXYBATE": "GABA-B受体、谷氨酸脱羧酶，中枢抑制",
    "OXYBATE": "GABA-B受体、谷氨酸脱羧酶，中枢抑制",
    "MAGNESIUM OXYBATE": "GABA-B受体、谷氨酸脱羧酶，中枢抑制",
    "CALCIUM OXYBATE": "GABA-B受体、谷氨酸脱羧酶，中枢抑制",
    "AMINOLEVULINIC ACID": "血红素合成通路，光动力治疗的光敏剂前体",
    "AMINOLEVULINIC ACID HYDROCHLORIDE": "血红素合成通路，光动力治疗的光敏剂前体",
    # 血液系统药物
    "TRIENTINE": "铜离子螯合剂，促进铜离子排泄",
    "TRIENTINE TETRAHYDROCHLORIDE": "铜离子螯合剂，促进铜离子排泄",
    "TRIENTINE HYDROCHLORIDE": "铜离子螯合剂，促进铜离子排泄",
    "AMIFOSTINE": "自由基清除剂，保护正常细胞免受放化疗损伤",
    "DEXTROSE": "葡萄糖转运体(GLUT)，细胞能量代谢底物",
    # 电解质/补充剂
    "SODIUM ACETATE": "碳酸氢根前体，酸碱平衡调节靶点",
    "POTASSIUM ACETATE": "钠-钾ATP酶，电解质平衡调节靶点",
    "CALCIUM ACETATE": "肠道磷结合靶点，钙敏感受体(CaSR)",
    "MAGNESIUM ACETATE": "钠-钾ATP酶、钙通道，电解质平衡调节",
    "ZINC ACETATE": "锌指蛋白、金属酶活性中心，肠道铜吸收抑制剂",
    "ALUMINUM ACETATE": "收敛靶点，蛋白质凝固，抗炎抗菌",
    # 消毒/防腐药
    "BENZALKONIUM CHLORIDE": "细菌细胞膜，破坏脂质双分子层",
    "CHLORHEXIDINE": "细菌细胞膜、细胞壁，破坏微生物完整性",
    # 其他
    "KRYPTON": "惰性气体，无特异性药理靶点，用于物理显像",
    "KRYPTON KR 81M": "惰性气体，无特异性药理靶点，用于物理显像",
    "SAPROPTERIN": "苯丙氨酸羟化酶辅酶，四氢生物蝶呤补充",
    "SAPROPTERIN DIHYDROCHLORIDE": "苯丙氨酸羟化酶辅酶，四氢生物蝶呤补充",
    "LISDEXAMFETAMINE": "去甲肾上腺素、多巴胺再摄取抑制剂，中枢神经兴奋",
    "LISDEXAMFETAMINE DIMESYLATE": "去甲肾上腺素、多巴胺再摄取抑制剂，中枢神经兴奋",
    "TEMOPORFIN": "光动力治疗光敏剂，产生活性氧损伤肿瘤细胞",
}

# --------------------------
# 通用工具函数
# --------------------------
def load_cache(cache_file):
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except:
            return {}
    return {}
def save_cache(cache_dict, cache_file):
    try:
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(cache_dict, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"⚠️ 缓存保存失败: {e}")
# --------------------------
# 核心API函数
# --------------------------
@lru_cache(maxsize=1000)
def get_smiles_by_id(ligand_id):
    if not ligand_id or ligand_id == "ZZZ":
        return None
    
    ligand_id = ligand_id.strip().upper()
    cache = load_cache(SMILES_CACHE_FILE)
    if ligand_id in cache:
        return cache[ligand_id]
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
        "Accept": "application/json"
    }
    
    max_retries = 2
    retry_delay = 1
    smiles_result = None
    # 1. PDBe v2 API
    for attempt in range(max_retries):
        try:
            url = f"https://www.ebi.ac.uk/pdbe/api/v2/compound/summary/{ligand_id}"
            r = requests.get(url, headers=headers, timeout=10)
            if r.status_code == 200:
                data = r.json()
                if ligand_id in data:
                    comp = data[ligand_id][0]
                    if "smiles" in comp:
                        for s in comp["smiles"]:
                            if s.get("name") == "canonical":
                                smiles_result = s.get("value")
                                break
                        if not smiles_result and comp["smiles"]:
                            smiles_result = comp["smiles"][0].get("value")
                        if smiles_result: break
            break
        except: time.sleep(retry_delay)
    # 2. RCSB PDB API
    if not smiles_result:
        for attempt in range(max_retries):
            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    if "rcsb_chem_comp_descriptor" in data:
                        desc = data["rcsb_chem_comp_descriptor"]
                        if "SMILES_stereo" in desc: smiles_result = desc["SMILES_stereo"]
                        elif "SMILES" in desc: smiles_result = desc["SMILES"]
                        if smiles_result: break
            except: time.sleep(retry_delay)
    # 3. PubChem XRef API
    if not smiles_result:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/PDB/{ligand_id}/cids/JSON"
            r = requests.get(url, headers=headers, timeout=10)
            if r.status_code == 200:
                cid = r.json()["IdentifierList"]["CID"][0]
                url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
                r2 = requests.get(url2, headers=headers, timeout=10)
                if r2.status_code == 200:
                    smiles_result = r2.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except: pass
    if smiles_result:
        cache[ligand_id] = smiles_result
        save_cache(cache, SMILES_CACHE_FILE)
        return smiles_result
    
    return None

# ====================== 【核心优化3：药物匹配逻辑新增靶点字段】======================
@lru_cache(maxsize=500)
def search_chembl_drugs(smiles, similarity_threshold):
    if not smiles: return [], "SMILES为空", 0
    
    cache_key = f"{smiles}_{similarity_threshold}"
    drug_cache = load_cache(DRUG_SEARCH_CACHE_FILE)
    if cache_key in drug_cache:
        return drug_cache[cache_key]["drugs"], drug_cache[cache_key]["msg"], drug_cache[cache_key]["total"]
    
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{safe_smiles}/{similarity_threshold}.json?limit=1000&expand=drug_indications"
        r = requests.get(url, headers=headers, timeout=30)
        
        if r.status_code == 200:
            mols = r.json().get('molecules', [])
            drugs = []
            
            for m in mols:
                if not (m.get('max_phase') and float(m.get('max_phase')) >= 4.0 and m.get('pref_name') and not m.get('withdrawn_flag')):
                    continue
                
                drug_name = m.get('pref_name', '').upper()
                
                # 智能适应症获取（优先级：手动字典 > API适应症 > 原型药匹配）
                indication_list = []
                if drug_name in MANUAL_DRUG_INDICATION:
                    indication_list.append(MANUAL_DRUG_INDICATION[drug_name])
                
                if not indication_list:
                    drug_indications = m.get('drug_indications', [])
                    for ind in drug_indications:
                        mesh_term = ind.get('mesh_heading', '').lower()
                        efo_term = ind.get('efo_term', '').lower()
                        for eng, cn in INDICATION_CN_MAP.items():
                            if eng in mesh_term or eng in efo_term:
                                indication_list.append(cn)
                                break
                
                if not indication_list:
                    base_name = drug_name.replace(' SULFATE', '').replace(' HYDROCHLORIDE', '').replace(' SODIUM', '').replace(' POTASSIUM', '').replace(' CALCIUM', '').strip()
                    if base_name in MANUAL_DRUG_INDICATION:
                        indication_list.append(MANUAL_DRUG_INDICATION[base_name])
                
                if not indication_list:
                    indication_list.append("暂无明确适应症")
                
                # 药物分类获取
                drug_class = "未分类"
                atc_list = m.get('atc_classifications', [])
                if atc_list:
                    atc_code = atc_list[0][:3]
                    if atc_code in ATC_FULL_MAP:
                        drug_class = ATC_FULL_MAP[atc_code]
                    else:
                        atc_first = atc_list[0][0].upper()
                        if atc_first in ATC_CLASS_MAP:
                            drug_class = ATC_CLASS_MAP[atc_first]
                
                # 【新增】官方靶点获取（优先级：手动字典 > 原型药匹配）
                drug_target = "暂无官方靶点"
                if drug_name in MANUAL_DRUG_TARGET:
                    drug_target = MANUAL_DRUG_TARGET[drug_name]
                else:
                    base_name = drug_name.replace(' SULFATE', '').replace(' HYDROCHLORIDE', '').replace(' SODIUM', '').replace(' POTASSIUM', '').replace(' CALCIUM', '').strip()
                    if base_name in MANUAL_DRUG_TARGET:
                        drug_target = MANUAL_DRUG_TARGET[base_name]
                
                # 【新增】官方靶点加入药物字典
                drugs.append({
                    "药物名称": m.get('pref_name'),
                    "ChEMBL ID": m.get('molecule_chembl_id'),
                    "相似度(%)": round(float(m.get('similarity', 0)), 2),
                    "分子量": round(float(m.get('molecule_properties', {}).get('full_mwt', 0)), 2),
                    "药物分类": drug_class,
                    "治疗适应症": "、".join(list(set(indication_list))),
                    "官方靶点": drug_target
                })
            
            drugs.sort(key=lambda x: x["相似度(%)"], reverse=True)
            
            drug_cache[cache_key] = {"drugs": drugs, "msg": "Success", "total": len(mols)}
            save_cache(drug_cache, DRUG_SEARCH_CACHE_FILE)
            
            return drugs, "Success", len(mols)
        else:
            return [], f"接口失败 ({r.status_code})", 0
    except Exception as e:
        return [], f"异常: {str(e)}", 0
# --------------------------
# 靶点加载函数
# --------------------------
def load_all_rna_targets():
    """加载所有RNA靶点"""
    print("🔍 正在加载所有RNA靶点...")
    
    try:
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        print("❌ 错误：未找到数据文件 'PDB_Dataset_Info_Full.xlsx'")
        exit(1)
    
    # RNA分类
    def categorize(desc):
        desc_lower = str(desc).lower()
        if 'riboswitch' in desc_lower: return "核糖开关 (Riboswitch)"
        elif 'aptamer' in desc_lower: return "适配体 (Aptamer)"
        elif any(word in desc_lower for word in ['quadruplex', 'g-4', 'g4']): return "G-四联体 (G-quadruplex)"
        elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']): return "核糖体 (rRNA)"
        elif 'ribozyme' in desc_lower: return "核酶 (Ribozyme)"
        elif any(word in desc_lower for word in ['ires', 'hairpin', 'stem-loop', 'pseudoknot']): return "特殊结构基元 (Special/Motifs)"
        else: return "其他 RNA (Others)"
    
    df['Category'] = df['Description (描述)'].apply(categorize)
    
    # 提取核心配体ID
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    
    # 筛选有配体的靶点
    valid_df = df[df['MainLigandID'] != "ZZZ"].copy()
    
    print(f"✅ 加载完成：共 {len(df)} 个PDB结构")
    print(f"   - 有核心配体的靶点: {len(valid_df)} 个")
    print(f"   - 无配体的靶点: {len(df)-len(valid_df)} 个（已跳过）")
    
    # 按类别统计
    category_stats = valid_df['Category'].value_counts()
    print(f"\n📊 靶点类别分布:")
    for cat, count in category_stats.items():
        print(f"   - {cat}: {count} 个")
    
    return valid_df
def load_antibacterial_targets():
    """加载抗菌RNA靶点（核糖开关和核糖体）"""
    print("🔍 正在加载抗菌RNA靶点...")
    
    try:
        df = pd.read_excel("PDB_Dataset_Info_Full.xlsx")
    except FileNotFoundError:
        print("❌ 错误：未找到数据文件 'PDB_Dataset_Info_Full.xlsx'")
        exit(1)
    
    # RNA分类
    def categorize(desc):
        desc_lower = str(desc).lower()
        if 'riboswitch' in desc_lower: return "核糖开关 (Riboswitch)"
        elif any(word in desc_lower for word in ['ribosomal', 'ribosome', 'rrna']): return "核糖体 (rRNA)"
        else: return "其他"
    
    df['Category'] = df['Description (描述)'].apply(categorize)
    
    # 提取核心配体ID
    def get_sort_id(ligand_text):
        if ligand_text == "No ligands" or str(ligand_text) == "nan": return "ZZZ"
        all_l = str(ligand_text).split(' | ')
        ions = ['MG', 'NA', 'K', 'CL', 'SO4', 'PO4', 'NCO', 'CD', 'ZN', 'CA', 'HG', 'FE', 'MN', 'CU', 'CO', 'BA', 'SR', 'RB', 'CS', 'LI', 'TL', 'BR', 'I', 'F']
        for l in all_l:
            lid = l.strip().split(' ')[0].upper()
            if lid not in ions: return lid
        return all_l[0].strip().split(' ')[0].upper()
        
    df['MainLigandID'] = df['Ligands (对应小分子)'].apply(get_sort_id)
    
    # 筛选抗菌靶点
    antibacterial_df = df[
        (df['Category'].isin(["核糖开关 (Riboswitch)", "核糖体 (rRNA)"])) & 
        (df['MainLigandID'] != "ZZZ")
    ].copy()
    
    # 优先筛选细菌相关靶点
    def is_bacterial(desc):
        desc_lower = str(desc).lower()
        return any(word in desc_lower for word in ['bacterial', 'bacteria', 'escherichia', 'coli', 'staphylococcus', 'aureus', 'streptococcus', 'mycobacterium', 'tuberculosis'])
    
    antibacterial_df['IsBacterial'] = antibacterial_df['Description (描述)'].apply(is_bacterial)
    antibacterial_df = antibacterial_df.sort_values(by=['IsBacterial', 'Category'], ascending=[False, True])
    
    print(f"✅ 筛选完成：共找到 {len(antibacterial_df)} 个抗菌RNA靶点")
    print(f"   - 核糖开关: {len(antibacterial_df[antibacterial_df['Category'] == '核糖开关 (Riboswitch)'])} 个")
    print(f"   - 核糖体: {len(antibacterial_df[antibacterial_df['Category'] == '核糖体 (rRNA)'])} 个")
    print(f"   - 明确细菌相关: {len(antibacterial_df[antibacterial_df['IsBacterial']])} 个")
    
    return antibacterial_df
# --------------------------
# 批量药物分析与报告生成
# --------------------------
def batch_analyze_targets(targets_df, similarity_threshold, output_dir, report_title, report_subtitle):
    """批量分析靶点并生成报告"""
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n🚀 开始批量药物重定位分析（相似度阈值: {similarity_threshold}%）...")
    
    all_results = []
    total_targets = len(targets_df)
    
    for idx, (_, row) in enumerate(targets_df.iterrows()):
        pdb_id = row['PDB ID']
        ligand_id = row['MainLigandID']
        category = row['Category']
        description = row['Description (描述)']
        
        # 每10个靶点打印一次进度
        if (idx + 1) % 10 == 0 or idx == 0:
            print(f"   正在处理 [{idx+1}/{total_targets}] {pdb_id} (配体: {ligand_id})")
        
        smiles = get_smiles_by_id(ligand_id)
        if not smiles:
            continue
        
        drugs, msg, total = search_chembl_drugs(smiles, similarity_threshold)
        
        for drug in drugs:
            drug["RNA类别"] = category
            drug["PDB ID"] = pdb_id
            drug["配体ID"] = ligand_id
            drug["RNA结构描述"] = description
            all_results.append(drug)
    
    print(f"\n✅ 批量分析完成！共匹配到 {len(all_results)} 条药物记录")
    
    if len(all_results) == 0:
        print("❌ 未匹配到任何药物")
        return
    
    # 优先级排序
    print("\n📊 正在进行药物优先级排序...")
    result_df = pd.DataFrame(all_results)
    
    # 计算优先级分数
    def calculate_priority(row):
        score = 0
        # 相似度权重：60%
        score += row["相似度(%)"] * 0.6
        # 药物分类权重：20%
        priority_classes = ["全身用抗菌药", "抗肿瘤药", "抗病毒药", "心血管系统药物"]
        if row["药物分类"] in priority_classes:
            score += 20
        # 特殊RNA类别权重：20%
        special_categories = ["核糖开关 (Riboswitch)", "核糖体 (rRNA)", "G-四联体 (G-quadruplex)"]
        if row["RNA类别"] in special_categories:
            score += 20
        return round(score, 2)
    
    result_df["优先级分数"] = result_df.apply(calculate_priority, axis=1)
    result_df = result_df.sort_values(by=["优先级分数", "相似度(%)"], ascending=False)
    
    # ====================== 【核心优化4：列顺序新增官方靶点】======================
    column_order = [
        "优先级分数", "RNA类别", "PDB ID", "配体ID", "药物名称", 
        "药物分类", "治疗适应症", "官方靶点", "相似度(%)", "分子量", "ChEMBL ID", "RNA结构描述"
    ]
    result_df = result_df[column_order]
    
    # 生成HTML报告
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 统计信息
    total_targets = len(targets_df)
    total_drugs = len(result_df)
    unique_drugs = result_df["药物名称"].nunique()
    # ====================== 【核心优化5：TOP20表格新增官方靶点列】======================
    top_20_drugs = result_df.head(20).to_html(index=False, classes="table table-striped table-hover")
    
    # 按RNA类别统计
    category_stats = result_df['RNA类别'].value_counts().to_dict()
    category_stats_html = "".join([f"<li><strong>{cat}</strong>: {count} 条药物记录</li>" for cat, count in category_stats.items()])
    
    # 按药物分类统计
    drug_class_stats = result_df['药物分类'].value_counts().head(10).to_dict()
    drug_class_stats_html = "".join([f"<li><strong>{cls}</strong>: {count} 条</li>" for cls, count in drug_class_stats.items()])
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="zh-CN">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{report_title}</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; }}
            .header {{ background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); color: white; padding: 2rem 0; margin-bottom: 2rem; }}
            .section {{ margin-bottom: 2rem; padding: 1.5rem; background: #f8f9fa; border-radius: 10px; }}
            .stat-card {{ background: white; padding: 1.5rem; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; }}
            .stat-number {{ font-size: 2rem; font-weight: bold; color: #11998e; }}
            table {{ font-size: 0.85rem; }}
            .footer {{ text-align: center; padding: 2rem; color: #6c757d; border-top: 1px solid #eee; margin-top: 3rem; }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="container">
                <h1 class="text-center">{report_title}</h1>
                <p class="text-center lead">{report_subtitle}</p>
                <p class="text-center">报告生成时间: {timestamp}</p>
            </div>
        </div>
        <div class="container">
            <div class="section">
                <h2>📊 分析摘要</h2>
                <div class="row">
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{total_targets}</div>
                            <div>分析靶点总数</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{total_drugs}</div>
                            <div>匹配药物总数</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{unique_drugs}</div>
                            <div>唯一药物数量</div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="stat-card">
                            <div class="stat-number">{len(category_stats)}</div>
                            <div>RNA类别数量</div>
                        </div>
                    </div>
                </div>
            </div>
            <div class="section">
                <h2>📈 数据分布统计</h2>
                <div class="row">
                    <div class="col-md-6">
                        <h5>RNA类别药物分布</h5>
                        <ul>{category_stats_html}</ul>
                    </div>
                    <div class="col-md-6">
                        <h5>药物分类TOP 10</h5>
                        <ul>{drug_class_stats_html}</ul>
                    </div>
                </div>
            </div>
            <div class="section">
                <h2>🏆 高优先级候选药物 TOP 20</h2>
                <p>基于相似度、药物分类和RNA类别综合排序 | 新增官方靶点字段</p>
                {top_20_drugs}
            </div>
            <div class="section">
                <h2>💡 研究建议</h2>
                <div class="row">
                    <div class="col-md-6">
                        <h5>🔬 抗菌药物方向</h5>
                        <ul>
                            <li>优先关注核糖开关和核糖体靶点</li>
                            <li>测试氨基糖苷类、大环内酯类药物的抗菌活性</li>
                            <li>针对超级细菌开发新型RNA靶向药物</li>
                        </ul>
                    </div>
                    <div class="col-md-6">
                        <h5>🎯 抗肿瘤药物方向</h5>
                        <ul>
                            <li>重点研究G-四联体靶点</li>
                            <li>筛选能够稳定G-四联体结构的药物</li>
                            <li>探索老药新用于肿瘤治疗的可能性</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>
        <div class="footer">
            <p>本报告由RNA结构精细化分类与AIDD药物重定位平台自动生成</p>
            <p>© 2024 RNA药物重定位自动化系统</p>
        </div>
    </body>
    </html>
    """
    
    report_path = os.path.join(output_dir, "drug_discovery_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    csv_path = os.path.join(output_dir, "candidate_drugs.csv")
    result_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    
    print(f"\n📄 报告生成完成！")
    print(f"   - 完整报告: {report_path}")
    print(f"   - 候选药物列表: {csv_path}")
    
    # 打印TOP 5
    print(f"\n📌 高优先级候选药物TOP 5:")
    for i in range(min(5, len(result_df))):
        row = result_df.iloc[i]
        print(f"   {i+1}. {row['药物名称']} (相似度: {row['相似度(%)']}%, RNA类别: {row['RNA类别']}, 靶点: {row['PDB ID']})")
    
    return result_df
# --------------------------
# 主函数
# --------------------------
def main():
    # 解析命令行参数
    mode = "all"  # 默认全RNA模式
    similarity_threshold = 40
    
    if len(sys.argv) > 1:
        mode = sys.argv[1].lower()
    if len(sys.argv) > 2:
        similarity_threshold = int(sys.argv[2])
    
    print("="*70)
    if mode == "antibacterial":
        print("🦠 全自动新型抗菌药物靶点发现与设计引擎")
        targets_df = load_antibacterial_targets()
        output_dir = ANTIBACTERIAL_OUTPUT_DIR
        report_title = "新型抗菌药物靶点发现与设计报告"
        report_subtitle = "基于核糖开关和核糖体RNA结构的AIDD药物重定位分析"
    else:
        print("🧬 全RNA靶点药物重定位自动化分析引擎")
        targets_df = load_all_rna_targets()
        output_dir = ALL_RNA_OUTPUT_DIR
        report_title = "全RNA靶点药物重定位分析报告"
        report_subtitle = "基于所有PDB RNA结构的系统性药物重定位分析"
    print("="*70)
    
    if len(targets_df) == 0:
        print("❌ 未找到任何有效靶点")
        return
    
    # 批量分析
    result_df = batch_analyze_targets(
        targets_df, 
        similarity_threshold, 
        output_dir, 
        report_title, 
        report_subtitle
    )
    
    print("\n" + "="*40)
    print("✅ 分析流程完成！")
    print("="*40)
    print(f"\n📂 所有结果已保存到 '{output_dir}' 目录")
if __name__ == "__main__":
    main()
