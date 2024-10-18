import numpy as np
import pyscf
from gau_pes import GauDriver

print("\nEnter PySCF!")
try:
    driver = GauDriver.from_stdio()

    mol = pyscf.gto.Mole()
    mol.spin = driver.multiplicity - 1
    mol.charge = driver.charge
    mol.basis = "sto-3g"
    mol.atom = driver.pyscf_atom()
    mol.build()

    mf = pyscf.scf.HF(mol)
    mf.chkfile = "h2o_pyscf.chk"
    mf.kernel()

    grad = None
    if driver.derivs in [1, 2]:
        grad = mf.nuc_grad_method().run().de

    force_constant = None
    if driver.derivs == 2:
        hess = mf.Hessian().kernel()
        atoms_mass = mol.atom_mass_list()
        atoms_mass_sqrt_inv = 1 / np.sqrt(atoms_mass)
        force_constant: np.ndarray = np.einsum(
            "ABij, A, B->ABij", hess, atoms_mass_sqrt_inv, atoms_mass_sqrt_inv)
        force_constant.swapaxes(1, 2)

    driver.write(
        energy=mf.energy_tot(),
        gradients=grad,
        force_constants=force_constant,
    )
except Exception as e:
    print(e)
print("Back to Gaussian!\n")
