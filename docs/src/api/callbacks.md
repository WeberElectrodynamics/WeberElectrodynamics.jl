# Callbacks

Callbacks are lightweight per-step hooks that fire around the algorithm's
core step. They can inspect or mutate integrator state (e.g. reflect a pair
that has passed through a collision singularity) without being fused into
the algorithm itself.

## Abstract interface

```@docs
HamiltonianCallback
apply_pre_step!
apply_post_step!
```

## Built-in callbacks

```@docs
CollisionBounce
```
