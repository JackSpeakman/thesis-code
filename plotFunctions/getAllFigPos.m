function pos = getNiceFigPos(layout,forPub)
% positions for figures (nice on my laptop)
%
% ------ Inputs -------
% layout        1-by-1          layout number
% forPub        1-by-1          for publication (binary)
%
% ------ Outputs -------
% pos           n_fig-by-4      position of figures


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
    otherwise
        pos = [];
end

end