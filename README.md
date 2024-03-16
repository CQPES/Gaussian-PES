# Gaussian-PES

Interface for invoking PES from Gaussian

## `BasePES`

`BasePES` is an abstract class that defines the protocol of a PES model. For some PES models with out analytical gradient, numerical gradient and force constant are also implemented.

## `GauDriver`

`GauDriver` is an I/O driver for parsing Gaussian internal input file and writing the results calculated on PES into the internal output file. Therefore, One can use `Gaussian` as an optimizer and does geometrical tasks, e.g. optimization, transition state, IRC, scan, and etc.

## Usage

Here is an example of defining a PES model and invoking it from Gaussian.

```python
# h2o_pes.py
from gau_pes import BasePES

_NUM_ATOMS = 3  # number of atoms of your system

class H2OPES(BasePES):
    def __init__(self) -> None:
        ...

    def calc_energy(
        self,
        coords: np.array,
    ) -> float
        """Calculate the potential energy of H2O system.

        Order of atoms: H H O
        """
        self._check_coords(_NUM_ATOMS, coords)

        # calculate energy
        energy = ...

        return energy


if __name__ == "__main__":
    from gau_pes import GauDriver

    driver = GauDriver()
    pes = H2OPES()

    driver.write(
        energy=pes.calc_energy(driver.coords),
        gradients=(pes.calc_gradients(driver.coords)
                   if driver.derivs in [1, 2] else None),
        force_constants=(pes.calc_force_constants(
            driver.coords) if driver.derivs == 2 else None),
    )
```

Then the input file be like

```
%nproc=1
%chk=H2O-opt-freq.chk
#p external="python3 -u /path/to/h2o_pes.py" opt=nomicro

H2O

0 1
 H     -0.00000002    0.00000000   -0.00000012
 H      0.92850569    0.00000000    1.19804603
 O     -0.00000008    0.00000000    0.95882669


--link1--
%nproc=1
%chk=H2O-opt-freq.chk
#p external="python3 -u /path/to/h2o_pes.py" freq geom=allcheck


```
