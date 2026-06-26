# Scene equation reference (manim MathTex, raw-string form)

Authoritative LaTeX for each scene's governing relation(s). Use with the house
style: `layout.equation(r"...")` (inline, accent colour) or
`layout.equation_card(r"...")` (the "hero" feature equation on a faint card).
Show feature equations reading-paced via `layout.show_equation(self, mob)`.
A "—" means the scene is visual; no new equation required (keep any it has).

## Preface
- **P00Hook** — —
- **P01Ellipse** — `r(\nu) = \dfrac{a(1-e^2)}{1 + e\cos\nu}` ; `q = a(1-e)` ; `Q = a(1+e)`
- **P02KeplerLaws** — `\dfrac{dA}{dt} = \tfrac12 r^2\dot\nu = \text{const}` ; `v^2 = GM\!\left(\dfrac{2}{r}-\dfrac{1}{a}\right)` ; `P^2 = a^3`
- **P03Elements** — `\varpi = \Omega + \omega` ; `M = E - e\sin E`
- **P04Geography** — `P \approx a^{3/2}\ \text{(yr, AU)}`
- **P05Clustering** — `\bar R = \left|\dfrac{1}{N}\sum_{k=1}^{N} e^{\,i\varpi_k}\right|`
- **P06Precession** — `\dot\varpi \approx \tfrac34\, n\, \dfrac{m_p}{M_\odot}\left(\dfrac{a_p}{a}\right)^{2}\dfrac{1}{(1-e^2)^2}`
- **P07Resonance** — `a_{\rm res} = a_p\left(\dfrac{p}{q}\right)^{2/3}` ; `\phi = (p+q)\lambda' - p\lambda - q\varpi`
- **P08Sculpting** — `\langle H\rangle = H(e,\,\Delta\varpi)` ; `\Delta\varpi \approx 180^\circ`
- **P09ReflectedLight** — `F \propto \dfrac{p\,R^2}{r^2\Delta^2} \approx \dfrac{1}{r^4}` ; `m = H + 5\log_{10}(r\,\Delta)`
- **P10Thermal** — `B_\nu(T) = \dfrac{2h\nu^3/c^2}{e^{h\nu/kT}-1}` ; `\lambda_{\rm peak} = \dfrac{b}{T}` ; `T_{\rm eq} = T_\odot\sqrt{\dfrac{R_\odot}{2d}}\,(1-A)^{1/4}`
- **P11Surveys** — `\text{detect} \iff m < m_{\rm lim}` ; `\Delta\theta = \mu\,\Delta t`
- **P12Indirect** — `\Delta\rho \propto \dfrac{GM_9}{d^3}` ; `\varepsilon_\odot \approx 6^\circ`
- **P13Map** — —

## Act I — evidence & clustering
- **Evidence2016** — `\bar R_\varpi` ; `P_{\rm joint} \approx 7\times10^{-5}`
- **Constraints2016** — `\bar R_{\Delta\varpi}` ; `f_{q>60\,\mathrm{AU}}`
- **Bias2017** — `N_{\rm intrinsic}(\varpi) = N_{\rm obs}(\varpi)/S(\varpi)`
- **OssosBias2017** — `N_{\rm det}(\varpi) \propto N_0\,S(\varpi)`
- **NapierCritique2021** — `p_{\rm Rayleigh} = e^{-N\bar R^2}`
- **ApsidalClustering2026** — `\bar R = \left|\tfrac1N\sum e^{\,i\varpi_k}\right|`
- **Clustering2019** — `\bar R,\quad p_{\rm Rayleigh}`
- **SheppardEtnos2016** — `\varpi_k \approx \langle\varpi\rangle`
- **NewDiscoveries2025** — `\Delta\varpi_k = \varpi_k - \langle\varpi\rangle`
- **Bp519Discovery2018** — `i \approx 54^\circ`
- **Clustering2025** — `\bar R(t) \sim \bar R_0\, e^{-t/\tau}`
- **PrimordialAlignment2024** — `\varpi(t) = \varpi_0 - \dot\varpi\,(t_0-t)`

## Act II — dynamics
- **SecularResonance2016** — `\dot\phi = 0,\quad \phi = \varpi - \varpi_9` ; `\ddot{\Delta\varpi} = -\omega_{\rm lib}^2\,\Delta\varpi`
- **Dynamics2017** — `\langle H\rangle(e,\Delta\varpi)` ; libration about `\Delta\varpi = 180^\circ`
- **SecularDynamics2018** — `(k,h) = e(\cos\Delta\varpi,\ \sin\Delta\varpi)` ; `e_{\rm forced}`
- **ChaoticDynamics2018** — `K = \dfrac{\delta a_{\rm res}}{\Delta a_{\rm sep}} \gtrsim 1`
- **InclinedTnos2016** — `\sqrt{1-e^2}\,\cos i \approx \text{const}`
- **DetachedInclinations2021** — `i_{\rm forced}(a)`
- **InclinationInstability2016** — `\gamma \sim n\,\dfrac{M_{\rm disk}}{M_\odot}`
- **Stability2021** — `q_{\rm crit}(a) = a_N\sqrt{\ln\!\left[\tfrac{576}{5}\tfrac{m_N}{M_\odot}\!\left(\tfrac{a}{a_N}\right)^{5/2}\right]}`
- **ResonancePrediction2016** — `a_{\rm res} = a_9\,(p/q)^{2/3}`
- **ResonanceHopping2017** — `a:\ \text{random walk among } a_{\rm res}`
- **Resonance2018** — `P_{\rm capture}(p{:}q)`
- **ResonanceHopping2020** — `\phi_{p:q} \ \text{librates}`
- **Commensurabilities2016** — `\dfrac{P_{\rm TNO}}{P_9} \approx \dfrac{p}{q}`
- **SecularOctupole2020** — `\epsilon_{\rm oct} = \dfrac{a}{a_9}\dfrac{e_9}{1-e_9^2}`
- **LowPerihelion2018** — `q_9 = a_9(1-e_9)`
- **Perturbation2025** — `\langle H\rangle = \langle H_0\rangle + \epsilon\langle H_1\rangle`
- **SelfgravDisk2019** — `\dot\varpi_{\rm disk} \propto \dfrac{M_{\rm disk}}{M_\odot}\,n`
- **OortSelfgrav2024** — `\dot\varpi_{\rm self} \propto \dfrac{M_{\rm IOC}}{M_\odot}\,n`
- **OortCloud2021** — `a \sim 2{,}000\text{–}20{,}000\ \mathrm{AU}`

## Act III — indirect signatures
- **CassiniRanging2016** — `\Delta\rho \propto GM_9\!\left[\dfrac{\mathbf r_9-\mathbf r_S}{|\mathbf r_9-\mathbf r_S|^3} - \dfrac{\mathbf r_9}{r_9^3}\right]`
- **HolmanPayne2016** — `m_{\max}(d) \propto d^{3}`
- **Obliquity2016** — `\varepsilon_\odot \approx 6^\circ` ; `\dot{\boldsymbol{s}} \propto \boldsymbol{s}\times\boldsymbol{\ell}_9`
- **ObliquityGomes2016** — `\varepsilon \to 6^\circ`
- **ObliquityLai2016** — `\varepsilon_\odot \propto \sin 2i_9`
- **Iorio2016** — `\dot\varpi_{\rm anom} \propto \dfrac{m_9}{d^3}`
- **Iorio2026** — `\dot\varpi_{\rm Saturn}^{\rm distended}`
- **UranusTilt2022** — `\varepsilon: 0^\circ \to 98^\circ`

## Act IV — the hunt (surveys)
- **WiseSearch2018** — `d_{\max}:\ m_{W1}(d_{\max}) = W1_{\rm lim}` ; `F_\nu = \pi B_\nu(T)(R/d)^2`
- **WiseCoadd2016** — `m_{\rm lim} \to m_0 + 1.25\log_{10}\!\sqrt{N}`
- **Ztf2021** — `f_{\rm excl} = \langle P_{\rm det}\rangle`
- **Des2022** — `f_{\rm DES}^{\rm unique} = \langle P_{\rm DES}(1-P_{\rm ZTF})\rangle`
- **PanStarrs2024** — `f_{\rm comb} = 1 - \prod_s (1 - P_s)`
- **TessYield2019** — `m_{\rm stack} = m_1 + 1.25\log_{10} N`
- **TessShiftstack2020** — `\Delta m = 1.25\log_{10} N`
- **Ps1Holman2025** — `\Delta m \approx 1\ \text{mag (shift–stack)}`
- **IrasAkari2025** — `\mu = \dfrac{\Delta\theta}{\Delta t}\ (23\ \mathrm{yr})`
- **IrasCandidate2022** — `F_{60} \Rightarrow d \sim \text{few}\times10^2\ \mathrm{AU}`
- **AkariRefutation2025** — `\text{orbit} \nRightarrow \text{consistency}`
- **SimonsForecast2025** — `d_{\max}:\ F_\nu^{\rm mm}(d_{\max}) = F_{\rm lim}`
- **ActMm2021** — `F_\nu \approx \dfrac{2k T \nu^2}{c^2}\dfrac{\pi R^2}{d^2}` (Rayleigh–Jeans)
- **LsstStrategy2023** — `f_{\rm disc} = f(\text{depth},\,N_{\rm link},\,\delta_{\min})`
- **ParallaxSearch2025** — `p = \dfrac{1\ \mathrm{AU}}{d}` ; `\Delta\theta_{\rm refl}`
- **Stacking2025** — `\mathrm{SNR} \propto \sqrt{N}`, trials penalty

## Act V — alternatives
- **MondEfe2023** — `a_0 \approx 1.2\times10^{-10}\ \mathrm{m\,s^{-2}}` ; external field `g_{\rm ext}`
- **PlanetY2025** — Laplace-plane warp `i_{\rm forced}(a)`
- **SirajOrbit2024** — posterior `P(m,a,e,i\mid \mathcal D)`
- **Orbit2021** — `\hat\theta_{\rm MCMC} = (6.2\,M_\oplus,\,380\,\mathrm{AU},\,0.3)`
- **Review2019** — `\approx 5\,M_\oplus,\ a\approx500\,\mathrm{AU},\ e\approx0.25`
- **CowanThermal2016** — `F_\nu = \pi B_\nu(T)(R/d)^2`
- **FortneyThermal2016** — `T_{\rm eff} = \max(T_{\rm eq},\,T_{\rm int})`
- **LinderEvolution2016** — `T_{\rm eff}(t) \downarrow,\quad L \propto R^2 T^4`
- **KuiperBelt2018** — `e = 1 - q/a`
- **RussellAlbedo2025** — `m = H + 5\log_{10}(r\,\Delta)` ; mini-Neptune `R = 2.0\text{–}2.6\,R_\oplus`, `p_V = 0.33\text{–}0.47`
- **StellarFlybys2025** — `\Gamma = n_\star\,\sigma\,v,\ \sigma = \pi q_\star^2\big(1+\tfrac{2GM}{q_\star v^2}\big)` ; `P \approx P_{\rm enc}\,f_{\rm geom}\,f_{\rm succ} \lesssim 5\%`
- **StellarCompanions2026** — `M_{\max}(d) = M_\star\,(d/d_\star)^3` ; halo `M(<1000\,\mathrm{AU}) < M_{\rm Pluto}`
- **AlphaSlope2026** — `F \propto d^{\alpha}`, `\alpha = -4` (reflected) vs `-2` (self-luminous)
