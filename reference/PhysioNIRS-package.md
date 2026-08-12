# PhysioNIRS: governed NIRS data and SNIRF input/output

`PhysioNIRS` maps continuous-wave SNIRF data into the
[`PhysioCore::PhysioExperiment()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/PhysioExperiment.html)
contract without changing measurement order, probe geometry, stimulus
tables, or the exact time base. Governed optical-density, modified
Beer-Lambert, motion-correction, short-separation nuisance,
quality/pruning, and live-neurofeedback operations retain their log,
distance, pathlength, geometry, filter, identity, and reference
conventions.

## Details

SNIRF metadata may contain identifying subject and acquisition fields.
Callers are responsible for de-identification before sharing a file.

## References

Shared Near Infrared Spectroscopy Format specification:
<https://fnirs.github.io/snirf/>

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioNIRS>

- <https://x-biosignal.r-universe.dev/PhysioNIRS>

- <https://x-biosignal.github.io/PhysioNIRS/>

- Report bugs at <https://github.com/x-biosignal/PhysioNIRS/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>

Authors:

- Yusuke Matsui <mail.to.matsui@gmail.com>
