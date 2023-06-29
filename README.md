# KEGG分析要点

# KEGG文件层级结构
- Pathway  
**1.** Pathway L1 (E.g Meatbolism)  
**2.** Pathway L2 (E.g Carbohydrate metabolism)  
**3.** pathway (E.g Glycolysis / Gluconeogenesis [PATH:ko00010])  
**4.** ko description  
  [Pathway层级解释](http://www.360doc.com/content/18/0313/11/45848444_736603164.shtml)  

- Module
Module 是KO的功能合集，每个module 代表1个功能单元  
[Module解释](https://www.jianshu.com/p/4c53fd4fc71c)

# KEGG数据预处理
- *KO的abundance*: 计算相对丰度（Cell Host Microbe. 2023 Feb 8;31(2):273-287.e5. doi: 10.1016/j.chom.2023.01.001.）  
- *ko的abundance*: 同属于一个ko的KO加合（Cell Host Microbe. 2022 Oct 12;30(10):1450-1463.e8. doi: 10.1016/j.chom.2022.09.004.）  

