# Particle in cell MIPT course.
Programs in chronological order:

1. nonrel-pusher.f95 &ndash; non-relativistic Boris particle pusher **for given fields** that is great at saving circles in constant magnetic fields while maybe not representing the phase correctly
1. pmove.py &ndash; program for drawing trajectory of single particle in 3d
2. rel-boris-pusher.f95 &ndash; relativistic Boris pusher
3. rel-vay-pusher.f95 &ndash; Vay pusher which is great at relativistic drift while being more expensive than Boris
4. NS-simulation.f95 &ndash; Vay pusher with some model neutron star fields
4. NS-sim.py &ndash; almost the same as pmove.py, but reading from other file
5. field-evol-2d.f95 &ndash; vacuum (no field and currents) field evolution program for 2d using Yee mesh, tested only with waves
5. field-evol-simple.f95 &ndash; program that draws profile of fields along x axis
6. pic-2d.f95 &ndash; full PIC program in 2d
