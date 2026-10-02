# Fe6+, Z=26, 20 electrons
nucmod = Nuclear.Model(26.0)
configs=[Configuration("[Ar] 3d^2"),Configuration("[Ar] 3d 4s"),Configuration("[Ar] 4s^2")]

grid = Radial.Grid(true)
ac=Atomic.Computation(Atomic.Computation(),nuclearModel=nucmod,grid=grid,configs=configs)
pac = perform(ac,output=true)
multiplet=pac["multiplet:"]
dl=Basics.displayLevels(stdout,[multiplet],N=1000)
