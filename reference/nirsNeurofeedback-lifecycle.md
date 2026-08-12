# Control and inspect live NIRS neurofeedback

`nirsNeurofeedbackStep()` processes at most `max_chunks` synchronously,
flushing committed target values after each chunk so the external update
timestamp remains the latest scope signal time.

## Usage

``` r
nirsNeurofeedbackStart(controller)

nirsNeurofeedbackStep(controller, max_chunks = 16L)

nirsNeurofeedbackState(controller)

nirsNeurofeedbackStop(controller)

nirsNeurofeedbackScope(controller)
```

## Arguments

- controller:

  A `NIRSNeurofeedback` runtime.

- max_chunks:

  Positive exact per-call work bound.

## Value

Lifecycle functions return `controller` invisibly. Step and state
functions return portable plain lists. `nirsNeurofeedbackScope()`
returns the governed `BiofeedbackScope` for a viewer.
