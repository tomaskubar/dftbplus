set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/fermihelper.F90
  ${curdir}/linearresponse.F90
  ${curdir}/numderivs2.F90
  ${curdir}/rotatedegen.F90
  ${curdir}/perturb.F90
  ${curdir}/finiteT.F90
  ${curdir}/perturbationx.F90
  ${curdir}/perturbationx_tk.F90
  ${curdir}/perturbationx_qmmm.F90
  ${curdir}/perturbationx_tk_qmmm.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
