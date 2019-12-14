# vasp_toolbox

The repository includes MATLAB codes & Linux codes used in VASP simulation (First principle calculation).

更新日志
2019.12.14
- 增加读取POSCAR文件（原胞中原子的个数）。
- 态密度文件DOSCAR现支持自旋极化/非极化(ISPIN  = 0/1)，自旋方向投影/轨道投影/不投影(LORBIT = 11/10/Default)的情况，相关数据读取到DOSCAR.dos和DOSCAR.projectedDOSforEachAtom元胞数组中。
已知缺陷：
- DOSCAR文件可以正确读入，但还没有制作不同情况下的绘制功能，需要绘制态密度可以选择DOSCAR.dos和DOSCAR.projectedDOSforEachAtom中的数据绘制。
- 绘制能带图时，内插点个数（KPOINTS文件第2行）不能随意更改（40、50经测试可用）。

2019.11.18
PlotBandAndDOS.m：
- 支持自旋极化情况下的能带图绘制。
- 支持自旋极化情况的自动识别（读取INCAR文件）。
已知缺陷：
- 自旋极化情况下，态密度文件读取错误，无法绘制正确的态密度分布。

MATLAB codes文件夹：

PlotBandAndDOS.m：
- 通过读取vasp计算的DOSCAR、EIGENVAL、KPOINTS、INCAR、POSCAR文件的信息，绘制态密度(DOS)和能带结构。
- 输入文件：能带计算得到的EIGENVAL、KPOINTS、INCAR、POSCAR和态密度计算得到的DOSCAR文件。将计算态密度得到的DOSCAR文件替换计算能带的DOSCAR文件后，选择整个文件夹作为输入文件也可。
- 输出结果：DOS图、能带结构图。
- 输入输出实例：PlotBandAndDOS_example文件夹中。
