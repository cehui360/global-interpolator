The data and code for the manuscript: HFSM: high-fidelity surface modeling for large-scale DEMs with terrain discontinuity preservation.
We share the publicly available airborne LiDAR benchmark dataset, provided by https://portal.opentopography.org. (../global-interpolator/data) for demonstration purposes.
A step-by-step instruction titled 'Instruction for model.docx' is provided (../global-interpolator/Instruction for model.docx), utilizing dataset to illustrate expected findings. 

1.Environments: 
-MATLAB R2022a

2.Description of data files:
In folder '../global-interpolator/data/':
-'Data_.txt': dataset.
Notice: Each row in the above tables represents x-coordinate, y-coordinate, and z-coordinate.

3.Description of codes: The code for model is stored in the '../global-interpolator/code' folder. 
-'VBCDA.m': Code for solution of the variational model objective function.
-'second_order_derivatives.m': Code for calculating the Dxv,Dyv,Dxxv,Dyyv,Dxyv.
-'Kronecker5.m': Code for calculating the Dx,Dy,Dxx,Dyy,Dxy.
-'saveFile.m': Code for saving the output file.
-'main.m': Code for starting the proposed method.





