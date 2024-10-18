#!/path/to/your/python/interpreter

import pyscf
from gau_pes import Orca6Driver

print("\nEnter PySCF!")
try:
    driver = Orca6Driver.from_stdio()

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
    if driver.do_gradient:
        grad = mf.nuc_grad_method().run().de

    driver.write(
        energy=mf.energy_tot(),
        gradients=grad
    )
except Exception as e:
    print(e)
print("Back to ORCA!\n")
