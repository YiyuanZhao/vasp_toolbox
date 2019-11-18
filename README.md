# vasp_toolbox

The repository includes MATLAB codes & Linux codes used in VASP simulation (First principle calculation).

更新日志
2019.11.18
PlotBandAndDOS.m：
- 支持自旋极化情况下的能带图绘制。
- 支持自旋极化情况的自动识别（读取INCAR文件）。
已知缺陷：
- 自旋极化情况下，态密度文件读取错误，无法绘制正确的态密度分布。

MATLAB codes文件夹：

PlotBandAndDOS.m：
- 通过读取vasp计算的DOSCAR、EIGENVAL、KPOINTS、INCAR文件的信息，绘制态密度(DOS)和能带结构。
- 输入文件：能带计算得到的EIGENVAL、KPOINTS、INCAR和态密度计算得到的DOSCAR文件。
- 输出结果：DOS图、能带结构图。
- 输入输出实例：PlotBandAndDOS_example文件夹中。
