function pos = getAllFigPosPub(caseStudy)
% positions for figures (standard so all plots look the same)

switch caseStudy
    case 'WO' 
        % Williams-Otto CSTR case study
        pos = [0,20,800,600;
            100,40,800,600;
            500,120,550,450;
            600,140,550,450;
            700,160,550,450];
    case 'WO2' 
        % Williams-Otto CSTR case study (2 constraint)
        pos = [0,20,550,450;
            100,40,550,450;
            200,60,550,450;
            500,120,550,450;
            600,140,550,450;
            700,160,550,450];
    otherwise
        pos = [];
end
        
end