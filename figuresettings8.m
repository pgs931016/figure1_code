function figuresettings8(filename, dpi, width, height)

    % 기본 설정
    alw = 0.5;    % Axes LineWidth
    fsz = 8;      % Fontsize
    lw=1;

    % 입력 인수 처리
    if nargin < 4
        height = 4; 
    end

    if nargin < 3
        width = 8; 
    end

    if nargin < 2
        dpi = 1200; 
    end

    if nargin < 1
        filename = ['unnamed_', datestr(datetime, 'yyyymmdd_HHMMSS')];
    end

    % Figure 설정
    set(0, 'DefaultLineLineWidth', 1);
    set(gcf, 'Units', 'centimeters');
    fig = gcf;
    fig.Position = [0, 0, width, height];
    set(gca, 'FontSize', fsz, 'LineWidth', alw, 'FontName', 'Times New Roman','LineWidth', lw);

    % 현재 축 가져오기
    ax = gca;

    % 축의 여백 정보 가져오기
    tightInset = ax.TightInset; 

    % 여백 설정
    leftMargin = tightInset(1);    % 왼쪽 여백
    bottomMargin = tightInset(2);  % 아래쪽 여백
    rightMargin = tightInset(3);   % 오른쪽 여백
    topMargin = tightInset(4);     

    % 축 내부 위치를 계산하여 균등하게 배치
    ax.Position = [leftMargin, bottomMargin, ...
                   1 - leftMargin - rightMargin, ...
                   1 - topMargin - bottomMargin];

    % 저장 설정
    set(gcf, 'InvertHardcopy', 'on');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [width, height]);
    set(gcf, 'PaperPosition', [0, 0, width, height]);

    % .fig 저장
    savefig([filename, '.fig']);

    % .tiff 저장
    print(gcf, [filename, '.tiff'], '-dtiff', ['-r', num2str(dpi)]);

end
% print('figure','-dpdf','-r2400');  
% 
% date = datetime('now','Format','yyyy-mm-dd');


 % idx = 1;
 %    while true
 %        filename = sprintf('%s_%s_%d', customname, date, idx);
 %        if ~isfile([filename '.pdf']) 
 %            break;
 %        end
 %        idx = idx + 1; 
 %    end
 % 
 %    % Save the figure as a PDF
 %    print(filename, '-dpdf', '-r2400');  


% if ispc % Use Windows ghostscript call
%   system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r2400 -oimprovedExample_eps.png improvedExample.eps');
% else % Use Unix/OSX ghostscript call
%   system('gs -o -q -sDEVICE=png256 -dEPSCrop -r2400 -oimprovedExample_eps.png improvedExample.eps');
% end
% end
