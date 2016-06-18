function out1 = plotsom_sahar(varargin)
  persistent INFO;
  if isempty(INFO), INFO = get_info; end
  if nargin == 0
    fig = nnplots.find_training_plot(mfilename);
    if nargout > 0
      out1 = fig;
    elseif ~isempty(fig)
      figure(fig);
    end
    return;
  end
  in1 = varargin{1};
  if ischar(in1)
    switch in1
      case 'info',
        out1 = INFO;
      case 'data_suitable'
        data = varargin{2};
        out1 = nnet.train.isNotParallelData(data);
      case 'suitable'
        [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
        [net,tr,signals] = deal(args{2:end});
        update_args = standard_args(net,tr,signals);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training_suitable'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals,param);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals);
        fig = nnplots.find_training_plot(mfilename);
        if isempty(fig)
          fig = figure('Visible','off','Tag',['TRAINING_' upper(mfilename)]);
          plotData = setup_figure(fig,INFO,true);
        else
          plotData = get(fig,'UserData');
        end
        set_busy(fig);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if isempty(unsuitable)
          set(0,'CurrentFigure',fig);
          plotData = update_plot(param,fig,plotData,update_args{:});
          update_training_title(fig,INFO,tr)
          nnplots.enable_plot(plotData);
        else
          nnplots.disable_plot(plotData,unsuitable);
        end
        fig = unset_busy(fig,plotData);
        if nargout > 0, out1 = fig; end
      case 'close_request'
        fig = nnplots.find_training_plot(mfilename);
        if ~isempty(fig),close_request(fig); end
      case 'check_param'
        out1 = ''; % TODO
      otherwise,
        try
          out1 = eval(['INFO.' in1]);
        catch me, nnerr.throw(['Unrecognized first argument: ''' in1 ''''])
        end
    end
  else
    [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
    update_args = standard_args(args{:});
    if ischar(update_args)
      nnerr.throw(update_args);
    end
    [plotData,fig] = setup_figure([],INFO,false);
    unsuitable = unsuitable_to_plot(param,update_args{:});
    if isempty(unsuitable)
      plotData = update_plot(param,fig,plotData,update_args{:});
      nnplots.enable_plot(plotData);
    else
      nnplots.disable_plot(plotData,unsuitable);
    end
    set(fig,'Visible','on');
    drawnow;
    if nargout > 0, out1 = fig; end
  end
end

function set_busy(fig)
  set(fig,'UserData','BUSY');
end

function close_request(fig)
  ud = get(fig,'UserData');
  if ischar(ud)
    set(fig,'UserData','CLOSE');
  else
    delete(fig);
  end
  drawnow;
end

function fig = unset_busy(fig,plotData)
  ud = get(fig,'UserData');
  if ischar(ud) && strcmp(ud,'CLOSE')
    delete(fig);
    fig = [];
  else
    set(fig,'UserData',plotData);
  end
  drawnow;
end

function tag = new_tag
  tagnum = 1;
  while true
    tag = [upper(mfilename) num2str(tagnum)];
    fig = nnplots.find_plot(tag);
    if isempty(fig), return; end
    tagnum = tagnum+1;
  end
end

function [plotData,fig] = setup_figure(fig,info,isTraining)
  PTFS = nnplots.title_font_size;
  if isempty(fig)
    fig = get(0,'CurrentFigure');
    if isempty(fig) || strcmp(get(fig,'NextPlot'),'new')
      if isTraining
        tag = ['TRAINING_' upper(mfilename)];
      else
        tag = new_tag;
      end
      fig = figure('Visible','off','Tag',tag);
      if isTraining
        set(fig,'CloseRequestFcn',[mfilename '(''close_request'')']);
      end
    else
      clf(fig);
      set(fig,'Tag','');
      set(fig,'Tag',new_tag);
    end
  end
  set(0,'CurrentFigure',fig);
  ws = warning('off','MATLAB:Figure:SetPosition');
  plotData = setup_plot(fig);
  warning(ws);
  if isTraining
    set(fig,'NextPlot','new');
    update_training_title(fig,info,[]);
  else
    set(fig,'NextPlot','replace');
    set(fig,'Name',[info.name ' (' mfilename ')']);
  end
  set(fig,'NumberTitle','off','ToolBar','none');
  plotData.CONTROL.text = uicontrol('Parent',fig,'Style','text',...
    'Units','normalized','Position',[0 0 1 1],'FontSize',PTFS,...
    'FontWeight','bold','ForegroundColor',[0.7 0 0]);
  set(fig,'UserData',plotData);
end

function update_training_title(fig,info,tr)
  if isempty(tr)
    epochs = '0';
    stop = '';
  else
    epochs = num2str(tr.num_epochs);
    if isempty(tr.stop)
      stop = '';
    else
      stop = [', ' tr.stop];
    end
  end
  set(fig,'Name',['Neural Network Training ' ...
    info.name ' (' mfilename '), Epoch ' epochs stop]);
end

%  BOILERPLATE_END  
%% =======================================================


function info = get_info
  info = nnfcnPlot(mfilename,'SOM Neighbor Distances combined with SOM Sample Hits',7.0,[]);
end

function args = training_args(net,tr,data)
  inputs = data.X;
  args = {net inputs};
end

function args = standard_args(varargin)
  [net,inputs] = deal(varargin{:});
  net = varargin{1};
  inputs = nntype.data('format',inputs);
  args = {net,inputs};
end

function plotData = setup_plot(fig)
  plotData.axis = subplot(1,1,1);
  plotData.numInputs = 0;
  plotData.numNeurons = 0;
  plotData.topologyFcn = '';
end


function fail = unsuitable_to_plot(param,net,input)
  if (net.numLayers < 1)
    fail = 'Network has no layers.';
  elseif (net.layers{1}.size == 0)
    fail = 'Layer has no neurons.';
  elseif isempty(net.layers{1}.distanceFcn)
    fail = 'Layer 1 does not have a distance function.';
  elseif isempty(net.layers{1}.topologyFcn)
    fail = 'Layer 1 does not have a topology function.';
  elseif ~strcmp(net.layers{1}.topologyFcn,'gridtop') ...
      && ~strcmp(net.layers{1}.topologyFcn,'hextop')
    fail = 'Only HEXTOP and GRIDTOP topology functions supported.';
  else
    fail = '';
  end
end

function plotData = update_plot(param,fig,plotData,net,inputs)

  inputs = inputs{1,1};
  numInputs = net.inputs{1}.processedSize;
  numNeurons = net.layers{1}.size;
  topologyFcn = net.layers{1}.topologyFcn;

  if strcmp(topologyFcn,'gridtop')  
    shapex = [-1 1 1 -1]*0.5;
    shapey = [1 1 -1 -1]*0.5;
    dx = 1;
    dy = 1;
    edgex = [-1 0 1 0]*0.5;
    edgey = [0 1 0 -1]*0.5;
  elseif strcmp(topologyFcn,'hextop')
    z = sqrt(0.75);
    shapex = [-1 0 1 1 0 -1]*0.5;
    shapey = [1 2 1 -1 -2 -1]*z/3;
    dx = 1;
    dy = sqrt(0.75);
    edgex = [-1 0 1 0]*0.5;
    edgey = [0 1 0 -1]*z/3;
  end
  shapex = shapex*0.3;
  shapey = shapey*0.3;

  pos = net.layers{1}.positions;
  dim = net.layers{1}.dimensions;
  numDimensions = length(dim);
  if (numDimensions == 1)
    dim1 = dim(1);
    dim2 = 1;
    pos = [pos; zeros(1,size(pos,2))];
  elseif (numDimensions > 2)
    pos = pos(1:2,:);
    dim1 = dim(1);
    dim2 = dim(2);
  else
    dim1 = dim(1);
    dim2 = dim(2);
  end

  
  if (plotData.numInputs ~= numInputs)  || any(plotData.dimensions ~= dim) || (plotData.numNeurons ~= numNeurons) ...
      || ~strcmp(plotData.topologyFcn,topologyFcn)
    set(fig,'NextPlot','replace');
    plotData.numInputs = numInputs;
    plotData.dimensions = dim;
    plotData.numNeurons = numNeurons;
    plotData.topologyFcn = topologyFcn;
    a = plotData.axis;
    set(a,...
      'DataAspectRatio',[1 1 1],...
      'Box','on',...
      'Color',[1 1 1],'FontSize',14) %background colour
    hold on

     plotData.p = zeros(1,numNeurons);
     plotData.t = zeros(1,numNeurons);
    % for i=1:numNeurons
    %   fill(pos(1,i)+shapex,pos(2,i)+shapey,[1 1 1], ...
    %     'EdgeColor',[0 0 0], ... %Edge colour of the neurons
    %     'FaceColor',[1 1 1]);   
    % end

    %Setup edges
    plotData.neighbors = sparse(tril(net.layers{1}.distances <= 1.001) - eye(numNeurons));
    plotData.numEdges = sum(sum(plotData.neighbors));
    plotData.patches = zeros(1,plotData.numEdges);
    plotData.text = zeros(1,plotData.numEdges);
    plotData.textt = zeros(1,plotData.numEdges);
    
    k = 1;
    for i=1:numNeurons
      for j=find(plotData.neighbors(i,:))
        pdiff = pos(:,j)-pos(:,i);
        angle = atan2(pdiff(2),pdiff(1));
        [ex,ey] = rotate_xy(edgex,edgey,angle);
        edgePos = (pos(:,i)+pos(:,j))*0.5;
        p1 = (2*pos(:,i) + 1*pos(:,j))./3;
        p2 = (1*pos(:,i) + 2*pos(:,j))./3;
        plotData.patches(k) = fill(edgePos(1)+ex,edgePos(2)+ey,[1 1 1],...
          'FaceColor',rand(1,3),...
          'EdgeColor','none');
        %plot([p1(1) p2(1)],[p1(2) p2(2)],'-','Color',[0 1 0]); %blue lines between neighbours
        k = k + 1;
      end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%writing the number of the region on dist_map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Setup neurons
   outputs = nncalc.y(net,{inputs});
  outputs = outputs{1};
  hits = sum(outputs,2);
 % norm_hits = sqrt(hits/max(hits));
 att=zeros(1,numNeurons);
 sim_t=sim(net, inputs);
for k1=1:numNeurons
    if (find(sim_t(k1,:)==1) > 0)
     at{k1}=find(sim_t(k1,:)==1);
    else
        at{k1}=[0];
    end
end

% for k1=1:numNeurons
%  att(1,k1)=at{1,k1};
% end



    for i=1:numNeurons
      fill(pos(1,i)+shapex,pos(2,i)+shapey,[1 1 1], ...
        'FaceColor',[188./255 143./255 188./255], ... %neuron colour is purple now
        'EdgeColor',[0 0 0]) %Edge colour of the neurons  
         if (hits(i) ~= 0) 
             plotData.t(i) = text(pos(1,i),pos(2,i),'', ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
        'FontWeight','bold',...
         'Color',[0 0 0], ...
         'FontSize',25);
         set(plotData.t(i),'String',num2str(at{1,i}));
         end
    end
    set(a,'XLim',[-1 (dim1-0.5)*dx + 1]);
    set(a,'YLim',[-1 (dim2-0.5)*dy + 0.5]);
    
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Done with writing the number of the region on dist_map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Finding weights between neurons and plot the colours between neurons.

  weights = net.IW{1,1};
  levels = zeros(1,plotData.numEdges);
  k = 1;
  for i=1:numNeurons
    for j=find(plotData.neighbors(i,:))
      levels(k) = sqrt(sum((weights(i,:)-weights(j,:)).^2));
      k = k + 1;
    end
  end
  mm = minmax(levels);
  levels = (levels-mm(1)) ./ (mm(2)-mm(1)); %adjust weights between 0 to 1
  if mm(1) == mm(2), levels = zeros(size(levels)) + 0.5; end

  k = 1;
  for i=1:numNeurons
    for j=find(plotData.neighbors(i,:))
      level = 1-levels(k);
%       red = min(level*2,1); % positive
%       green = max(level*2-1,0); % very positive/negative
      %c = [red green 0];   %colours between neurons
      %c = min(1,nngui.red*2*(1-level));
      c = [1 1 1]*(level);
     
      set(plotData.patches(k),'FaceColor',c,'EdgeColor',[0 0 0]);
%         pdiff = pos(:,j)-pos(:,i);
%         angle = atan2(pdiff(2),pdiff(1));
%         [ex,ey] = rotate_xy(edgex,edgey,angle);
%         edgePos = (pos(:,i)+pos(:,j))*0.5;
%         p1 = (2*pos(:,i) + 1*pos(:,j))./3;
%         p2 = (1*pos(:,i) + 2*pos(:,j))./3;
%         text1=int2str(level*100);
%          text(p1(1), p1(2),text1, ...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%         'FontWeight','bold',...
%          'Color',[1 0 0], ...
%          'FontSize',10);
        % set(plotData.text(k),'String',num2str(level*100));
      k = k + 1;
    end
  end
end

function [x2,y2] = rotate_xy(x1,y1,angle)
  [a,r] = cart2pol(x1,y1);
  a = a + angle;
  [x2,y2] = pol2cart(a,r);
end

