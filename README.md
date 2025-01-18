# Fractional-DMD

*FO-DMD* reconstructs the dynamic with multiple modes. The dynamic is defined as follows:
```math
f(x,t)=cos(x)\text{E}_{0.5,1} ((1+1i)t^{0.5})+\cos(2x+10)\text{E}_{0.5,1}((-2+1i)t^{0.5}),
```
where $x\in [-10,10]$ and $t\in[0,10]$

Run `FODMD.m` will generate the time response at $x = −10$, spatiotemporal dynamics of the real model and reconstructed models, and the spatial modes captured by *FO-DMD*

