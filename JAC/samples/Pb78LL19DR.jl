# Compute the low-lying resonances of the DR spectrum of Pb.78+ (1s^2 2s^2 --> 1s·2 2s 2p ^3p_1) ions # An electron is captured into shells with n >= 19 
setDefaults("print summary: open", "Pb78LL19DR.jac")
setDefaults("unit: energy", "eV") 
setDefaults("unit: rate", "1/s") 
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox=15.0) 
asfSettings = AsfSettings(AsfSettings(), scField=Basics.DFSField(1.0)) # not used 
iConfigs = [Configuration("1s^2 2s^2")] 
eConfigs = [Configuration("1s^2 2s 2p")] 
fConfigs = [Configuration("1s^2 2s^2 2p")] 
addShells = Basics.generateShellList(19, 19, [l for l=0:11]) 
toShells = Basics.generateShellList(2, 5, [0, 1, 2, 3, 4]) 
finalConfigs = Basics.generateConfigurations(fConfigs, [Shell("2s"), Shell("2p")] , toShells) 
intermConfigs = Basics.generateConfigurationsWithAdditionalElectron(eConfigs, addShells) 
intermConfigs = Basics.merge(intermConfigs,finalConfigs) 
drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings(), multipoles = [E1],printBefore=true, gauges = [JAC.UseCoulomb, JAC.UseBabushkin], calcOnlyPassages=true, 
corrections = DielectronicRecombination.AbstractCorrections[DielectronicRecombination.HydrogenicCorrections(18,missing,missing), DielectronicRecombination.ResonanceWindowCorrection(0., 2.0)])
comp = Atomic. Computation(Atomic. Computation(), name="Low-lying DR resonances for Pb^78+", grid=grid, nuclearModel=Nuclear.Model(82.), initialConfigs = iConfigs, intermediateConfigs = intermConfigs, finalConfigs = 
finalConfigs, processSettings = drSettings)
results = perform(comp; output=true)
setDefaults("print summary: close", "")


