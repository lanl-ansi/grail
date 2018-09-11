# grail

Gas Reliability Analysis Integrated Library: algorithms for natural gas pipeline optimization, optimal control, and simulation

# Instructions:

1. Add the IPOPT MEX binary to the root folder, available at https://www.coin-    or.org/download/binary/Ipopt/ or included with the Opti toolbox at     https://www.inverseproblem.co.nz/OPTI/

2.    grail.exe can be compiled using the MATLAB compiler, and can be distributed for use with     MATLAB runtime, which is available here: https://www.mathworks.com/products/compiler/matlab-    runtime.html

3.    Specify the folder containing the model and parameters in model_folder.txt (e.g. model30t) and add it to the root. Running grail will read the input in that folder and put output in the same folder.


    These publications can be cited as references:

Zlotnik, A., Chertkov, M. and Backhaus, S.,  Optimal control of transient flow in natural gas networks. In 2015 IEEE 54th Annual Conference on Decision and Control (CDC),  (pp. 4563-4570). IEEE (2015)

Zlotnik, A., Roald, L., Backhaus, S., Chertkov, M. and Andersson, G.,  Coordinated scheduling for interdependent electric power and natural gas infrastructures. IEEE Transactions on Power Systems, 32(1), pp.600-610. (2017)

Sundar, K. and Zlotnik, A., State and Parameter Estimation for Natural Gas Pipeline Networks using Transient State Data. arXiv preprint arXiv:1803.07156. (2018)