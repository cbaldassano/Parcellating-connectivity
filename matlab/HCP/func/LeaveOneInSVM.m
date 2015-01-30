function [LOI_libsvm, LOI_libsvm_models] = LeaveOneInSVM(parcelData, patientData_Y, parcellation)

% parcelData = cell(length(patientData_X),1);
% numParcels = max(parcellation);
% for s = 1:length(patientData_X)
%     parcelData{s} = zeros(size(patientData_X{s},1), numParcels);
%     for i = 1:numParcels
%         parcelData{s}(:,i) = mean(patientData_X{s}(:, parcellation == i),2);
%     end
%     parcelData{s} = parcelData{s} ./ repmat(std(parcelData{s}), size(parcelData{s},1), 1);
% end

LOI_libsvm = zeros(length(parcelData),1);
LOI_libsvm_models = cell(length(parcelData),1);
parfor i = 1:length(parcelData)
    disp(['Leaving in patient ' num2str(i)]);

    % Concatenate all training data except patient i
    Xsv_train = parcelData{i};
    Ysv_train = patientData_Y{i};
    Xsv_test = vertcat( parcelData{1:length(parcelData) ~= i} );
    Ysv_test = vertcat( patientData_Y{1:length(patientData_Y) ~= i} );

    LOI_libsvm_models{i} = svmtrain(Ysv_train, Xsv_train, '-t 0');

    [~, acc, ~] = svmpredict(Ysv_test, Xsv_test, LOI_libsvm_models{i});
    LOI_libsvm(i) = acc(1);
 end
        

end