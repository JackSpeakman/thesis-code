function pos = getAllFigPos(caseStudy)
% positions for figures (nice on my laptop)

switch caseStudy
    case 'WO'
        % Williams-Otto CSTR case study
        pos = [0,320,640,320;
            640,320,640,320;
            0,20,430,230;
            430,20,430,230;
            860,20,430,230];
    case 'WO2'
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