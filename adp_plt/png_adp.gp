#!/usr/bin/gnuplot

# plot
set terminal png enhance
set output "F.png"
set xlabel "rho / unknown"
#set xlabel "{/Symbol r} / a.u."
set ylabel "Embedding function, F(rho) / eV"
#set ylabel "Embedding Function, F({/Symbol r}) / eV"
plot "F.plt" u 1:2 w l t "F(rho)"
#plot "F.plt" u 1:2 w l t "F({/Symbol r})"
#
set output "rho.png"
set xlabel "r / Angstrom"
set ylabel "Density function, rho(r) / a.u."
#set ylabel "Density Function, {/Symbol r}(r) / a.u."
plot "rho.plt" u 1:2 w l t "rho(r)"
#plot "rho.plt" u 1:2 w l t "{/Symbol r}(r)"
#
set output "rphi.png"
set xlabel "r / Angstrom"
#set ylabel "Effective charge function, Z(r) / (Hartree*Bohr-radii)^0.5"
set ylabel "Effective Charge Function, r*phi(r) / eV*Angstrom"
plot "rphi.plt" u 1:2 w l t "rphi(r)"
#
set output "u.png"
set xlabel "r / Angstrom"
set ylabel "dipole function / a.u."
#set ylabel "Density Function, {/Symbol r}(r) / a.u."
plot "u.plt" u 1:2 w l t ""
#
set output "w.png"
set xlabel "r / Angstrom"
set ylabel "Quadrupole function / a.u."
#set ylabel "Density Function, {/Symbol r}(r) / a.u."
plot "w.plt" u 1:2 w l t " "
#
#set output "diff_energy.png"
#set xlabel "Step, N"
#set ylabel "Differential Energy, E - E(ref) / eV"
#plot "energy.dat" u 1:2 w l t "E-E(ref)"
