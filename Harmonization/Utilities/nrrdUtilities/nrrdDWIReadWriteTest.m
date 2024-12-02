
%% Function to test the reading/writing of NRRD files
%% The purpose of this file is to test if nrrd files can
%% be read/written without loss of information
function [ errorStatus ] = nrrdDWIReadWriteTest( )
    testDataNames = {  'Synthetic-Test-data_01.nhdr', 'Synthetic-Test-data_02.nhdr', 'Synthetic-Test-data_03.nhdr' };
    
    for testDataIndex=1:size(testDataNames,2)
        testBaseName = char(testDataNames(testDataIndex));
        testFn = fullfile(fileparts(pwd),'TestSuite',testBaseName);

        % Check that the file exists
        assert(exist(testFn, 'file') == 2, 'File does not exist');

        errorStatus = 0;
        % Test current teem based Nrrd File readers/writers
        errorStatus = errorStatus + testNrrdReaderWriters( testFn );

        % Test new ITK based readers/writers
        errorStatus = errorStatus + testITKReaderWriters( testFn );
    end
end


function [ errorValue ] = testNrrdReaderWriters ( filename1 )
    sprintf('** testNrrdReaderWriters: %s **',filename1)
    firstPassRaw = nrrdLoadWithMetadata(filename1);
    secondPassFn = fullfile(tempdir(),'firstout.nhdr');
    nrrdSaveWithMetadata(secondPassFn,firstPassRaw);
    secondPassRaw = nrrdLoadWithMetadata(secondPassFn);

    errorValue = compareDWIdata( firstPassRaw, secondPassRaw );
    errorStatus = ~ ( ( isfinite(errorValue)) &&  (errorValue < sqrt(eps)) );
    if errorStatus 
        sprintf('ERROR: Teem Test failed structures are different!')
    else
        sprintf('SUCCESS!: Teem Test passed structures are the same!')
    end
    sprintf('%f,%d\n%s\n%s',errorValue,errorStatus,filename1,secondPassFn)
end


function [ errorValue ] = testITKReaderWriters ( filename1 )
    sprintf('** testITKReaderWriters: %s **',filename1)
    firstPassRaw = itkLoadWithMetadata(filename1);
    secondPassFn = fullfile(tempdir(),'firstITKout.nhdr');
    nrrdSaveWithMetadata(secondPassFn,firstPassRaw);
    secondPassITKFn = fullfile(tempdir(),'firstITKout-2.nhdr');
    itkSaveWithMetadata(secondPassITKFn,firstPassRaw);
    secondPassRaw = nrrdLoadWithMetadata(secondPassFn);

    errorValue = compareDWIdata( firstPassRaw, secondPassRaw );
    errorStatus = ~ ( ( isfinite(errorValue)) &&  (errorValue < sqrt(eps)) );
    if errorStatus 
        sprintf('ERROR: Teem Test failed structures are different!')
    else
        sprintf('SUCCESS!: Teem Test passed structures are the same!')
    end
    sprintf('%f,%d\n%s\n%s',errorValue,errorStatus,filename1,secondPassFn)
end

function [ errorValue ] = compareDWIdata( dwi1, dwi2 )
    errorValue = 0.0;
    errorValue = errorValue + compareImages( dwi1, dwi2 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'bvalue',0.0 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'gradientdirections',0.0001 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'measurementframe', 0.0 );
end

function [ errorValue ]  = compareImages ( im1, im2 )
    errorValue = 0.0;
    errorValue = errorValue + abs(double(im1.space)-double(im2.space));
    errorValue = errorValue + compareFloatArrays(im1, im2, 'data', 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'spacedirections', 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'spaceorigin', 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'centerings', 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'kinds',0.0 );  % Check the enumerations are same

    errorValue = errorValue + compareStringArrays( im1, im2, 'spaceunits');
    errorValue = errorValue + compareStringArrays( im1, im2, 'spacedefinition');
end

function [ error ] = compareFloatArrays( fa, sa, fieldName, tolerance )
  %% Compare first and second array to determine how different they are
  if isfield(fa, fieldName) && isfield(sa, fieldName)
      array1 = double(fa.(fieldName));
      array2 = double(sa.(fieldName));
      array1size = size(array1)
      array2size = size(array2)
      denom = norm(array2(:))^2;
      if ~( isfinite(denom) && denom > eps )
          denom = 1;
      end
      error = norm(array1(:)-array2(:))^2/denom;
  else
      sprintf('ERROR: MISSING FIELD: %s',fieldName);
      error = 10000000;
  end
  if ( ~isfinite(error) ) ||  (error > tolerance)
      sprintf('ERROR: Differences found for: %s, (%f > %f)', fieldName, error, tolerance)
  end
end

function [ error ] = compareStringArrays( fa, sa, fieldName )
  %% Compare first and second array to determine how different they are
  if isfield(fa, fieldName) && isfield(sa, fieldName)
      sarray1 = fa.(fieldName);
      sarray2 = sa.(fieldName);
      error = sum(1.-strcmp(sarray1,sarray2));
  else
      sprintf('ERROR: MISSING FIELD: %s',fieldName);
      error = 10000000;
  end
  if error > 0
      sprintf('ERROR: Differences found for: %s', fieldName)
  end
end
