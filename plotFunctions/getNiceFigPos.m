function pos = getNiceFigPos(layout,forPub)
% getNiceFigPos gets the positions for figures (nice on my laptop)
%
% ------ Inputs -------
% layout        1-by-1          layout number
% forPub        1-by-1          for publication (binary)
%
% ------ Outputs -------
% pos           n_fig-by-4      position of figures

if forPub
    % for publication
    switch layout
        case 1
            % Williams-Otto CSTR case study
            pos = [0,20,550,450;
                100,40,550,450;
                500,120,550,450;
                600,140,550,450;
                700,160,550,450];
        case 2
            % Williams-Otto CSTR case study (2 constraint)
            pos = [0,20,550,450;
                100,40,550,450;
                200,60,550,450;
                500,120,550,450;
                600,140,550,450;
                700,160,550,450];
        case 3
            % Distillation column case study
            pos = [0,20,550,450;
                300,60,550,450;
                600,120,550,450];
    end
else
    % for screen
    switch layout
        case 1
            % Williams-Otto CSTR case study
            pos = [0,320,640,320;
                640,320,640,320;
                0,20,430,230;
                430,20,430,230;
                860,20,430,230];
        case 2
            % Williams-Otto CSTR case study (2 constraint)
            pos = [0,350,430,300;
                430,350,430,300;
                860,350,430,300;
                0,20,430,300;
                430,20,430,300;
                860,20,430,300];
        case 3
            % Distillation column case study
            pos = [0,320,640,350;
                640,320,640,350;
                320,20,640,350];
    end
end

end