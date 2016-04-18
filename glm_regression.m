function [s_bold_t_series_phys, yDelta, adjRsq ] = glm_regression(s_bold_t_series, design_matrix, order, dicom_header, BW, HW_flag, retroicor_flag)

%Calculate yDelta and Radj^2 

yDelta = zeros(dicom_header.ydim, dicom_header.xdim, dicom_header.tdim);

s_bold_t_series_phys = zeros(size(yDelta));

adjRsq = zeros(dicom_header.ydim, dicom_header.xdim);

design_matrix(isnan(design_matrix)) = 0; 

t_i=(1:1:dicom_header.tdim)';

T_design_matrix = array2table(design_matrix,'VariableNames',{'ac1', 'bc1', 'ac2', 'bc2', 'ac3', 'bc3', 'ac4','bc4','ac5','bc5', 'ar1', 'br1', 'ar2', 'br2', 'ar3', 'br3','ar4','br4','ar5','br5', 'x111','x112','x113','x114','x211','x212','x213','x214','x121','x122','x123','x124','x221','x222','x223','x224'});

if retroicor_flag == 1
    % Modify design_matrix and T_design_matrix based on specified order
    if order.X == 0
        if order.C == 0
            T_design_matrix_mod = T_design_matrix(:,[11:(10+2*order.R)]); % just R terms
        elseif order.R == 0
            T_design_matrix_mod = T_design_matrix(:,[1:(2*order.C)]); % just C terms
        else
            T_design_matrix_mod = T_design_matrix(:,[1:(2*order.C),11:(10+2*order.R)]); % just CR terms
        end
    else
        if order.C == 0
            T_design_matrix_mod = T_design_matrix(:,[11:(10+2*order.R),21:(20+4*order.X)]); % just RX terms
        elseif order.R == 0
            T_design_matrix_mod = T_design_matrix(:,[1:(2*order.C),21:(20+4*order.X)]); % just CX terms
        else
             T_design_matrix_mod = T_design_matrix(:,[1:(2*order.C),11:(10+2*order.R),21:(20+4*order.X)]); % CRX terms
        end
    end

    design_matrix_mod = table2array(T_design_matrix_mod);
end

h1=waitbar(0,'Performing linear regression...');

for i=1:dicom_header.ydim
    
    for j=1:dicom_header.xdim
        
        if BW(i,j) == 1
                
            y=squeeze(s_bold_t_series(i,j,:));

            test_signal = y - mean(y);

            T_sig = table(test_signal,'VariableNames',{'BOLD_Signal'});
            
            if HW_flag == 1
            
                p1 = polyfit(t_i,test_signal,1);
            
                HW = p1(2) + p1(1)*t_i;
                
            elseif HW_flag == 2
                
                p2 = polyfit(t_i,test_signal,2);
            
                HW = p2(3) + p2(2)*t_i + p2(1)*t_i.^2;
            
            end
            
            T_HW = array2table(HW, 'VariableNames', {'HW2'});
            
            if retroicor_flag == 1
            
                T_regressors = [T_HW T_design_matrix_mod];

                ext_design_matrix = [ones(dicom_header.tdim,1) HW design_matrix_mod];       
                
            else
                
                T_regressors = T_HW;
                
                ext_design_matrix = [ones(dicom_header.tdim,1) HW];
                
            end
                
            T_glm = [T_regressors T_sig];

            mdl = fitlm(T_glm,'linear');

            coef_tbl = mdl.Coefficients(:,1);

            coefvals = table2array(coef_tbl)';

            coefvals = repmat(coefvals, dicom_header.tdim ,1);

            yDelta(i,j,:) = sum(coefvals.*ext_design_matrix,2);
            
            adjRsq(i,j) = mdl.Rsquared.Adjusted;
            
            s_bold_t_series_phys(i,j,:) = y - squeeze(yDelta(i,j,:));
        
        else
            s_bold_t_series_phys(i,j,:) = s_bold_t_series(i,j,:);
            
            adjRsq(i,j) = NaN;
        end
    end

    waitbar(double(i)/double(dicom_header.ydim), h1);

end

delete(h1);


