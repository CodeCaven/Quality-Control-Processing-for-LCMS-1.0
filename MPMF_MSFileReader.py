from pymsfilereader import MSFileReader
import numpy as np
import matplotlib.pyplot as plt

filename = "HelaiRT1uL_190722200502"
filename2 = "QC_Metabolomics_190707160342"
rawfile = MSFileReader("Z:\\qc_automation\\fusion\\instrument_data\\" + filename + ".raw")
# Z:\\Metabolomics\\QC_runs\\C1_Clayton
# Z:\\qc_automation\\fusion\\instrument_data\\

#print('Version', rawfile.Version())
#print('GetFileName', rawfile.GetFileName())
#print('GetCreatorID', rawfile.GetCreatorID())
#print('GetVersionNumber', rawfile.GetVersionNumber())
#print('GetCreationDate', rawfile.GetCreationDate())
#print('IsError', rawfile.IsError())
#print('IsNewFile', rawfile.IsNewFile())
#print('IsThereMSData', rawfile.IsThereMSData())
#print('HasExpMethod', rawfile.HasExpMethod())
#print('InAcquisition', rawfile.InAcquisition())
#print('GetErrorCode', rawfile.GetErrorCode())
#print('GetErrorMessage', rawfile.GetErrorMessage())
#print('GetWarningMessage', rawfile.GetWarningMessage())
#print('RefreshViewOfFile', rawfile.RefreshViewOfFile())
#print('GetNumberOfControllers', rawfile.GetNumberOfControllers())

#print("GetNumberOfControllersOfType('No device')", rawfile.GetNumberOfControllersOfType('No device'))
#print("GetNumberOfControllersOfType('MS')", rawfile.GetNumberOfControllersOfType('MS'))
#print("GetNumberOfControllersOfType('Analog')", rawfile.GetNumberOfControllersOfType('Analog'))
#print("GetNumberOfControllersOfType('A/D card')", rawfile.GetNumberOfControllersOfType('A/D card'))
#print("GetNumberOfControllersOfType('PDA')", rawfile.GetNumberOfControllersOfType('PDA'))
#print("GetNumberOfControllersOfType('UV')", rawfile.GetNumberOfControllersOfType('UV'))
##print("GetControllerType('MS')", rawfile.GetControllerType('MS'))


#print('GetCurrentController()', rawfile.GetCurrentController())
#print('GetExpectedRunTime()', rawfile.GetExpectedRunTime())
#print('GetMaxIntegratedIntensity()', rawfile.GetMaxIntegratedIntensity())
#print('GetMaxIntensity()', rawfile.GetMaxIntensity())
#print('GetInletID()', rawfile.GetInletID())
#print('GetErrorFlag()', rawfile.GetErrorFlag())
#print('GetFlags()', rawfile.GetFlags())
#print('GetAcquisitionFileName()', rawfile.GetAcquisitionFileName())
#print('GetOperator()', rawfile.GetOperator())
#print('GetComment1()', rawfile.GetComment1())
#print('GetComment2()', rawfile.GetComment2())
#print('GetFilters()', rawfile.GetFilters())
#print('GetMassTolerance()', rawfile.GetMassTolerance())

#print('GetMassResolution', rawfile.GetMassResolution())
#print('GetNumTrailerExtra', rawfile.GetNumTrailerExtra())
#print('GetLowMass', rawfile.GetLowMass())
#print('GetHighMass', rawfile.GetHighMass())
#print('GetStartTime', rawfile.GetStartTime())
#print('GetEndTime', rawfile.GetEndTime())
#print('GetNumSpectra', rawfile.GetNumSpectra())
#print('GetFirstSpectrumNumber', rawfile.GetFirstSpectrumNumber())
#print('GetLastSpectrumNumber', rawfile.GetLastSpectrumNumber())
#print('GetAcquisitionDate', rawfile.GetAcquisitionDate())
#print('GetUniqueCompoundNames', rawfile.GetUniqueCompoundNames())

#print('############################################## INSTRUMENT BEGIN')
#print('GetInstrumentDescription', rawfile.GetInstrumentDescription())
#print('GetInstrumentID', rawfile.GetInstrumentID())
#print('GetInstSerialNumber', rawfile.GetInstSerialNumber())
#print('GetInstName', rawfile.GetInstName())
#print('GetInstModel', rawfile.GetInstModel())
#print('GetInstSoftwareVersion', rawfile.GetInstSoftwareVersion())
#print('GetInstHardwareVersion', rawfile.GetInstHardwareVersion())
#print('GetInstFlags', rawfile.GetInstFlags())
#print('GetInstNumChannelLabels', rawfile.GetInstNumChannelLabels())
# #print( 'GetInstChannelLabel(0)', rawfile.GetInstChannelLabel(0) )
#print('IsQExactive', rawfile.IsQExactive())  # Not implemented in MSFileReader 3.0.29.0
#print('############################################## INSTRUMENT END')

#print('GetChroData', rawfile.GetChroData(startTime=rawfile.StartTime,
                                            #endTime=rawfile.EndTime,
                                            #massRange1="{}-{}".format(rawfile.LowMass, rawfile.HighMass),
                                            #scanFilter="Full ms "))
# get and plot UV controllers
uv_number = rawfile.GetNumberOfControllersOfType('UV')
for i in range(1, uv_number+1):
    rawfile.SetCurrentController(4,i) # 4 is UV enumerated number type, controllers start at 1
    print('Plotting controller..', rawfile.GetCurrentController())

    uv_data = rawfile.GetChroData(startTime=rawfile.StartTime, endTime=rawfile.EndTime)
    x = np.array(uv_data[0][0])
    y = np.array(uv_data[0][1])
    print(len(x))

    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=1)
    ax.set(xlabel='Time (minutes)', ylabel='Pressure', title='UV Controller ' + str(i))
    fig.savefig("UV_controller_" + str(i) + "_" + filename + ".png")
    plt.show()

rawfile.Close()