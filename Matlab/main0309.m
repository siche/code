%% * main program 
% process video
% deal the 
% preparetion


cd('D:/20180309');
dbstop if error;

% set detection parameter
gap = 150;
i = 0;
sig_num=1;
conv_factor=ones(2);
conv_num=3;
sig_size=1;

% min cluster distance
distance=16; 

% stop point used in test 
output_num=80000;

% exclude two signal in one block
b2_flag=0;
l2_flag=0;
r2_flag=0;


% flag is the number of pixels that has signal
left_flag=[];
right_flag=[];
left=[];right=[];

% detection boundary
lrow1=15;lrow2=85;
lcol1=45;lcol2=116;

rrow1=21;rrow2=89;
rcol1=131;rcol2=210;

bg_range=min(lrow1,lcol1);


% initialize result 
n1=0;n2=0;n3=0;
n12=0;n22=0;
l=0;r=0;
ll=0;rr=0;bb=0;
l1=[];r1=[];
l2=[];r2=[];b2=[];
% prepartion for output image
% if there is no dir result, will make dir
% else will remove present dir and make new dir

 if (~exist('result','dir'))
    mkdir ('result/1/');
    mkdir ('result/0/');
    mkdir ('result/l2/');
    mkdir ('result/r2/');
    mkdir ('result/b2/');
 else
    rmdir result s;
    mkdir ('result/1/');
    mkdir ('result/0/');
    mkdir ('result/l2/');
    mkdir ('result/r2/');
    mkdir ('result/b2/');
end
%}

% num of video 
n=21;

for k=1:n
% break the output
%
if i>=output_num 
    break
end
%}

disp([num2str(k),'/',num2str(n)]);

video_file =['video/delay9ns_full_pump_100w_X',num2str(k),'.avi'];
video = VideoReader(video_file);
% hei = video.Height;
% len = video.Width;
% sig = zeros(hei,len);
% process frmaes in video
while hasFrame(video)
    i=i+1;
    p=readFrame(video);
    p=rgb2gray(p);
    background = p(1:bg_range,1:bg_range);
    b=mean(mean(background));
    threshold=gap+b;

    sig=double(p)>threshold;
    sig_conv=conv2(double(sig),conv_factor);

    % sig1=conv2(sig,conv_factor);
    
    % check if there is signal in every signal block
    sig1=sig(lrow1:lrow2,lcol1:lcol2);
    sig1_conv=sig_conv(lrow1:lrow2,lcol1:lcol2)>=conv_num;
    % ceil round towards positive infinity
    sig1=ceil(sum(sum(sig1_conv))/sig_size);
    left_flag(i)=sig1;

    sig2=sig(rrow1:rrow2,rcol1:rcol2);
    sig2_conv=sig_conv(rrow1:rrow2,rcol1:rcol2)>=conv_num;
    sig2=ceil(sum(sum(sig2_conv))/sig_size); 
    right_flag(i)=sig2;
    
    % output image
    %{
    if i<=output_num
        if(sig1>=sig_num || sig2>=sig_num)
            output_image(p,i,1);
       % else
       %     output_image(p,i,0);
        end
    else 
        break
    end
    %}

    % calculate n1 n2 n3...
    
    
    % L Block
    if (sig1>=sig_num && sig2<sig_num)
        
        % judge if there is two signals
        %if sig1>2*sig_num
        l2_flag=0;
        loc=find(sig1_conv);
        point_num=length(loc);
        if (point_num>2 && sig1>=2*sig_num)
            [rows,cols]=size(sig1_conv);

            row=mod(loc,rows);
            col=floor(loc/rows)+1;
            row(row==0)=rows;

            location=[row,col];
	        [idx,c]=kmeans(location,2);
            dis=c(1,:)-c(2,:);
            dis=dis*dis';
            
            if dis>distance
                n12=n12+1;
                ll=ll+1;
                l2(ll)=i;
                l2_flag=1;

                output_image(p,i,2);
            end
        end
        
        if(~l2_flag)
            l=l+1;
            n1=n1+1;
            l1(l)=i;
            output_image(p,i,1);
        end
        
    end

    
 % R Block
    if(sig1<sig_num && sig2>=sig_num)
        % output_image(p,i,1);
        %if sig2>2*sig_num
        r2_flag=0;
        loc=find(sig2_conv);
        point_num=length(loc);
        if (point_num>2 && sig2>=2*sig_num)
            [rows,cols]=size(sig2_conv);
            row=mod(loc,rows);
            row(row==0)=rows;
            col=floor(loc/rows)+1;

            location=[row,col];
	        [idx,c]=kmeans(location,2);
            dis=c(1,:)-c(2,:);
            dis=dis*dis';
	 
            if dis>distance
                n22=n22+1;
                rr=rr+1;
                r2(rr)=i;
                r2_flag=1;
                output_image(p,i,3);
            end
        end
        
        if(~r2_flag)
             n2=n2+1;
             r=r+1;
	         r1(r)=i;
             output_image(p,i,1)
        end
    end

  
 % Both Block
    b2_flag=(l2_flag | r2_flag);
    if (sig1>=sig_num && sig2>=sig_num && ~b2_flag)
        n3=n3+1;
	    bb=bb+1;
	    b2(bb)=i;
        output_image(p,i,4);
    end
end
end
%% * calculate n1,n2,n3

%{
n1=0;n2=0;n3=0;
n12=0;n22=0;

l=0;r=0;b=0;
ll=0;rr=0;
l1=[];r1=[];
l2=[];r2=[];b2=[];
% n1:left 
% n2:right
% n3:both left and right
% n4:left and right at the same location
% n12 left has more than two signals
% n22 right has more than two signals 

for j=1:i
    sig1=left_flag(j);
    sig2=right_flag(j);

    if (sig1>=sig_num && sig2<sig_num)
        n1=n1+1;
	l=l+1;
	l1(l)=j;
    end

    if(sig1<sig_num && sig2>=sig_num)
        n2=n2+1;
	r=r+1;
	r1(r)=j;
    end

    % set the threshold as 1
    if sig1>=2*sig_num        
        n12=n12+1;
	ll=ll+1;
	l2(ll)=j;
    end

    if sig2>=2*sig_num
        n22=n22+1;
	rr=rr+1;
	r2(rr)=j;
    end

    % set threshold as 2
 
    if (sig1>=sig_num && sig2>=sig_num)
        n3=n3+1;
	b=b+1;
	b2(b)=j;
    end
end
%}

message=['sig_num=',num2str(sig_num),...
	 '\n conv_num=',num2str(conv_num),...
	 '\n sig_size=',num2str(sig_size),...
     '\n frames=',num2str(i)];
% mail2me('sichecmm@sjtu.edu.cn','Result of 20180206',...
%    'Running result','result.mat');

save result.mat *
mfile=[mfilename,'.m'];
zip('result',{mfile,'result.mat','result'});
mail2me('sichecmm@sjtu.edu.cn','result of 20180206',message,'result.zip');


function output_image(img,frame_num,option)
% output signaled image to specified path
    if option==2
        filename = ['result/l2/',num2str(frame_num),'.jpg'];
    elseif option==3
        filename =['result/r2/',num2str(frame_num),'.jpg'];
    elseif option==4
        filename = ['result/b2/',num2str(frame_num),'.jpg'];
    else
        filename = ['result/',num2str(option),'/',num2str(frame_num),'.jpg'];
    end
    imwrite(img, filename);
end

