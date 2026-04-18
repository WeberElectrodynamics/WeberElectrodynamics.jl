# Wolfram Language scripts

Sandbox for running symbolic/numeric Wolfram Language computations alongside the
Julia package. Runs on the free [Wolfram Engine](https://www.wolfram.com/engine/)
via `wolframscript`.

## Run

```bash
wolframscript -file hello.wls
```

Or, since `hello.wls` has a shebang:

```bash
chmod +x hello.wls
./hello.wls
```

## First-time setup

After installing Wolfram Engine, activate the free Developer License once:

```bash
wolframscript -activate
```

Sign in with a free [Wolfram ID](https://www.wolfram.com/id/).

## Files

- `hello.wls` — smoke test: symbolic integration, high-precision numerics,
  linear algebra, and a simple ODE via `DSolve`.
