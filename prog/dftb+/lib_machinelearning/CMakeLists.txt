set(curdir "lib_machinelearning")

set(sources-fpp
  ${curdir}/MachineLearning.F90
  ${curdir}/MachineLearningNeuralNet.F90
  ${curdir}/MachineLearningSymmetryFunctions.F90)


set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
