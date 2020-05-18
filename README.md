### ["Can Global Uncertainty Promote International Trade?"](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/bvw_june2019.pdf)

Code repository for ["Can Global Uncertainty Promote International Trade?"](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/bvw_june2019.pdf) by By [Isaac Baley](http://www.isaacbaley.com/), [Laura Veldkamp](https://www.waugheconomics.com/), and [Michael Waugh](https://www.waugheconomics.com/)

Journal of International Economics, May 2020

---

**PRELIMINARIES**

- Some routines use functions of the compecon library by Miranda and Fackler.
   [http://www4.ncsu.edu/~pfackler/compecon/toolbox.html](http://www4.ncsu.edu/~pfackler/compecon/toolbox.html)

- See Appendix C for details on computation and parameterization.  

- Color code:  BLUE (MATLAB code), GREEN (MATLAB matrix or Excel), VIOLET (EPS figures)

---
**Figure 1: Trade Policy Uncertainty and Exports**

The folder ``Figure-1`` contains:

- [``TPUandExports.xlsx``](Figure-1/TPUandExports.xlsx) Excel file with two data series.  
1.	US exports from NIPA (quarterly frequency).

2.	US Trade Policy Uncertainty (monthly frequency). Categorical Economic Policy Uncertainty Index for the U.S. constructed by Baker, Bloom, and Davis and downloaded from http://www.policyuncertainty.com. Variables are normalized to 100 for 2014Q1.

- [``Figure1.m``](Figure-1/Figure1.m)

 - Reads Excel data from [``TPUandExports.xlsx``](Figure-1/TPUandExports.xlsx)

 - Plots fig1_tpuexports.eps

---

**Figure 2: Export Coordination under Perfect and Imperfect Information**

**Figure 3: Uncertainty Increases Terms of Trade Volatility and Average**

The folder [``Figures-2-and-3``](Figures-2-and-3/) contains:

- [``Simu_twoprecisions.m``](Figures-2-and-3/Simu_twoprecisions.m):
 Solves the model with CES preferences under perfect (zero signal noise) and imperfect information (limit to infinite signal noise) and simulates for T=20 periods.

  - Calls [``PI.m``](Figures-2-and-3/PI.m) to guess initial coefficients from solution of perfect information case.

  - Calls [``modelsolve.m``](Figures-2-and-3/modelsolve.m) to solve for imperfect information policy.

  - Calls [``simulation.m``](Figures-2-and-3/simulation.m) to generate series of outcomes (exports and terms of trade).

- Saves the results from the simulations in [``results_twoprecisions.mat``](Figures-2-and-3/results_twoprecisions.mat)

- [``Figure23.m``](Figures-2-and-3/Figure23.m)

  - Reads [``results_twoprecisions.mat``](Figures-2-and-3/)

  - Plots ``fig2_simutradelevel.eps``

  - Plots ``fig3_simutermstrade.eps``

---
**Figure 4: Effect of Uncertainty in Trade Depends on the Elasticity of Substitution**

The folder [``Figure 4``](Figures-4/) contains:

- [``Solve_twoelasticities.m``](Figures-4/Solve_twoelasticities.m)  Solves the model for two elasticities of substitution and various levels of signal noise.
  - Calls [``PI.m``](Figures-4/PI.m)  to guess initial coefficients from solution of perfect information case.  

  - Calls [``modelsolve.m``](Figures-4/modelsolve.m)  to solve for imperfect information policy.

  - Saves the results in [``results_twoelasticities.mat``](Figures-4/results_twoelasticities.mat)

- [``Figure4.m``](Figures-4/Figure4.m)

  - Reads results_twoelasticities.mat

  - Plots fig4_elasticity.eps

---

**Figure 5: Higher Uncertainty Brings Average Foreign Beliefs Closer to the Prior**

The folder [``Figure-5``](Figure-5/) contains:

- [``Figure5.m``](Figure-5/Figure5.m) Computes average foreign beliefs for low and high realizations of domestic endowment, for various levels of signal noise.

  - Plots ``fig5_priorbeliefs.eps``

---

**Figure 6: State-dependent Responses to Uncertainty**

The folder [``Figure-6``](Figure-6/) contains:

- [``Solve_statedep.m``](Figure-6/Solve_statedep.m) c Solves the model for low elasticities of substitution and various levels of signal noise.

  - Calls [``PI.m``](Figure-6/PI.m) c to guess initial coefficients from solution of perfect information case.

  - Calls [``modelsolve.m``](Figure-6/modelsolve.m) c to solve for imperfect information policy.

  - Saves the results in results_statedep.mat

- [``Figure6.m``](Figure-6/Figure6.m) c

  - Reads results_statedep.mat

  - Plots fig6_statedep.eps

---
**Figure 7: Uncertainty Decreases Trade Coordination and Increases Utility Correlation**

The folder [``Figure-7``](Figure-7/) contains:

- [``Solve_risksharing.m``](Figure-7/Solve_risksharing.m) Solves the model for a low elasticity of substitution and various levels of signal noise.

   - Calls [``PI.m``](Figure-7/PI.m) to guess initial coefficients from solution of perfect information case.  

   - Calls [``modelsolve.m``](Figure-7/modelsolve.m) to solve for imperfect information policy.

   - Calls [``simulation.m``](Figure-7/simulation.m) to compute trade and utility correlation across T=100,000 draws.

   - Saves the results in [``results_statedep.mat``](Figure-7/results_statedep.mat)

- [``Figure7.m``](Figure-7/Figure7.m)

  - Reads ``results_risksharing.mat``

  - Plots ``fig7_risksharing.eps``

---

**Figure 8: Completing the Market Reduces Exports**

The folder [``Figure-8``](Figure-8/) contains:

- [``Solve_financial.m``](Figure-8/Solve_financial.m)  Solves the model with two types of agents: a fraction with perfect information (complete markets) and a fraction with imperfect information. Considers low elasticity of substitution.

  - Calls [``PI.m``](Figure-8/PI.m)  to guess initial coefficients from solution of perfect information case.  

  - Calls [``modelsolve_financial.m``](Figure-8/modelsolve_financial.m)  to compute equilibrium policies.

  - Calls [``simulation_financial.m``](Figure-8/simulation_financial.m)  to simulate T=1,000 draws.  

  - Saves the results in ``results_financial.mat``

- [``Figure8.m``](Figure-8/Figure8.m)

  - Reads ``results_financial.mat``

  - Plots ``fig8_financial.eps``

---

**Figure 9 (Appendix): Comparative statistics for risk aversion**

The folder [``Figure-9``](Figure-9/) contains:

- [``Solve_twoelasticities_riskaversion.m``](Figure-9/Solve_twoelasticities_riskaversion.m)  Solves the model three times with different levels of risk aversion. For each level of risk aversion, solves the model for low and high elasticities of substitution and various levels of signal noise.

  - Calls [``PI.m``](Figure-9/PI.m)  to guess initial coefficients from solution of perfect information case.

  - Calls [``modelsolve.m``](Figure-9/modelsolve.m)  to solve for imperfect information policy.

  - For sigma = 1-theta (benchmark), saves the results in [``results_ra_bench.mat``](Figure-9/)

  - For sigma = 0 (neutral), saves the results in [``results_ra_neutral.mat``](Figure-9/)

  - For sigma = 1.5 (averse), saves the results in [``results_ra_averse.mat``](Figure-9/)

- [``Figure9.m``](Figure-9/Figure9.m)

  - Reads the three files: ``results_ra_bench.mat``, ``results_ra_neutral.mat``, and  ``results_ra_averse.mat``

  - Plots ``fig9_riskaversion.eps``
