clc;
clear ;
close ;

%% 导入文件
load('eeg_9.mat');

%% 常规参数
fs=128;%采样率
T=10;%取窗的时间长度
epoch=T*fs;
overlap=0.8*epoch;%正常交叠窗取75%；不要超过0.8


%% 特殊参数
r=0.15;%阈值,fast_ApEn和sampen

%% 伪参数
d=10;
D=d*fs;% number of FFT points.frequency resolution=0.5
warning off;
%% 读取数据
for No=1%（1~58,60岁包含,59~100,70岁包含）%length(TianjingData)
    EEGch1=a.pEEG1(1,:);
    %     bis=TianjingData(No).BIS;
    %     name=TianjingData(No).name;
    %     time=double(TianjingData(No).Time);
    
    length_Raweeg1=length(EEGch1);
    N=floor((length_Raweeg1-epoch)/(epoch-overlap))+1;%生成指标长度
    
    %分配内存
    %     BIS=resample(bis,N,length(bis));%对bis降采样
    %     Time=resample(time,N,length(time));%对bis降采样
    pe=zeros(N,1);
    apen=zeros(N,1);
    sampen=zeros(N,1);
    SFS=zeros(N,1);
    Delta=zeros(N,1);
    Theta=zeros(N,1);
    Alpha=zeros(N,1);
    Beta=zeros(N,1);
    Gamma=zeros(N,1);
    RE=zeros(N,1);
    SE=zeros(N,1);
    MPF=zeros(N,1);
    SEF95=zeros(N,1);
    plzc=zeros(N,1);
    dfa=zeros(N,1);
    %     Physio=zeros(N,4);
    %% 生成指标
    for p=1:N
        
        dataSegments=EEGch1((p-1)*(epoch-overlap)+1:(p-1)*(epoch-overlap)+epoch);
        
        pe(p,1)=perm(dataSegments,6,1);%根据经验，fs=10,epoch=10*fs，m=6,tao=1
        %         apen(p,1)=fast_ApEn(dataSegments,r);%【2EEG entropy……(P4)】 from liang:N=1000,r=0.1~0.25,m=2~3,fast_ApEn equal to ApEn(when m=2)
        %         sampen(p,1)=fast_SampEn(dataSegments,2,r);%邵帅:dataSegments_N=10*100,m=2,r=0.15;【2EEG entropy……(P5)】:基于ApEn的改进，取值等同于ApEn
        [SFS(p,1),Delta(p,1),Theta(p,1),Alpha(p,1),Beta(p,1),Gamma(p,1)] = sfsAndPower(dataSegments,epoch,T);%参数被T确定
        %         [SEF95(p,1),MPF(p,1)] = sef95AndMPF(dataSegments,epoch,T);%参数被T确定
        %         [RE(p,1),SE(p,1)] = SEandRE(dataSegments,T);%参数被T确定
        plzc(p,1)=PLZC(dataSegments,4,1);%邵帅:dataSegments_N=1000（10*100）,m=4,tao=1
        %         dfa(p,1)=DFA_a(dataSegments);%DFA_a调用DFA最终算出a；a=0.5无自相关，a为(0,0.5)正相关,(0.5,1)负相关，a=1信号为噪声
        %         Physio(p,:)=Phy(No,:);%生理参数
        %         BIS(p,1)=bis(p);%重采样后的BIS
        %         Time(p,1)=Time(p);%重采样后的BIS
        p
    end
    sampleSet{No,1}=[pe,SFS,Delta,Theta,Alpha,Beta,Gamma,plzc];
    save sampleSet_pEEG1 sampleSet%每生成1个就保存一个
%     clear pe apen sampen SFS Delta Theta Alpha Beta Gamma SEF95 MPF RE SE plzc dfa
end
