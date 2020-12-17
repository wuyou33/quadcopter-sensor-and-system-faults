classdef quadVisualizationManager
    %quadVisualizationManager To viauslize multiple quadcopter
    %   Detailed explanation goes here
    
    properties
        quadData %= quadVisualization
        quadNum
        quadLine
        quadColor = {'red','green','yellow','magenta','cyan','black'}
        quadPathColor = {'black','cyan','blue','magenta','yellow','green','red'}
        quadPathLineType = {'--','-.','-','-*'}
    end
    
    methods
        function obj = quadVisualizationManager(quadNum)
            %UNTITLED9 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                obj.quadData = quadVisualization;
                obj.quadNum = 1;
            else
                obj.quadNum = quadNum;
                for ii = 1:quadNum
                    obj.quadData(ii) = quadVisualization;
                end
            end
        end
        
         function obj = visualize(obj,xpath,ypath,psipath,mass,dt)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % Plot wind velocity
            xmin = -5;xmax = 15;
            ymin = -5;ymax = 15;
            [X,Y] = meshgrid(xmin:3:xmax,ymin:3:ymax);
            U = 1+0.5*rand(size(X));
            V = -1+0.5*rand(size(Y));
            f1 = figure(1);clf
            f1.Position = [-1500 150 1280 800];
            h = subplot(3,2,[1 2 3 4]);
            h.XLim = [-5 15];
            h.YLim = [-5 15];
            h.XLabel.String = 'X [m]';
            h.YLabel.String = 'Y [m]';
            hold on
            g = quiver(X,Y,U,V);
            
            %Create the quadcopter frame
            ii = 1;
            lineWidth = 1.5;
            for jj = 1:obj.quadNum
                first_color = obj.quadColor{jj};
            obj.quadData(jj) = obj.quadData(jj).update(xpath(ii,jj),ypath(ii,jj),psipath(ii,jj));

            %Create line objects
            obj.quadLine(jj).path = line(xpath(ii,jj),ypath(ii,jj),'Color',obj.quadPathColor{jj},'linestyle',obj.quadPathLineType{jj},'linewidth',lineWidth);
            obj.quadLine(jj).x(1) = line(obj.quadData(jj).propeller(:,1,1),obj.quadData(jj).propeller(:,2,1),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).x(2) = line(obj.quadData(jj).propeller(:,1,2),obj.quadData(jj).propeller(:,2,2),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).x(3) = line(obj.quadData(jj).propeller(:,1,3),obj.quadData(jj).propeller(:,2,3),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).x(4) = line(obj.quadData(jj).propeller(:,1,4),obj.quadData(jj).propeller(:,2,4),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).x(5) = line(obj.quadData(jj).frame(:,1,1),obj.quadData(jj).frame(:,2,1),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).x(6) = line(obj.quadData(jj).frame(:,1,2),obj.quadData(jj).frame(:,2,2),'Color',first_color,'linewidth',lineWidth);
            obj.quadLine(jj).orientation = line(obj.quadData(jj).orient(:,1),obj.quadData(jj).orient(:,2),'Color',first_color,'linewidth',lineWidth);
            end
            h2 = subplot(3,2,5);
            h2.XLim = [0 120];
            h2.YLim = [-1 1];
            h2.XLabel.String = 'Time [s]';
            h2.YLabel.String = 'Yaw angle [degree]';
            hold on
            yawangle = line([0 dt],[0 psipath(1)],'Color',first_color,'linewidth',lineWidth);
            h3 = subplot(3,2,6);
            h3.XLim = [0 120];
            h3.YLim = [0.4 1.5];
            h3.XLabel.String = 'Time [s]';
            h3.YLabel.String = 'm [Kg]';
            massLine = line([0 dt],[mass(1) mass(1)],'Color',first_color,'linewidth',lineWidth);
            
            % Plot in real time
            Npath = length(xpath);
            t1 = tic;
            
            for ii = 1:Npath
                for jj = 1:obj.quadNum
                obj.quadData(jj) = obj.quadData(jj).update(xpath(ii,jj),ypath(ii,jj),psipath(ii,jj));
                
                obj.quadLine(jj).x(1).XData = squeeze(obj.quadData(jj).propeller(:,1,1));
                obj.quadLine(jj).x(1).YData = squeeze(obj.quadData(jj).propeller(:,2,1));
                obj.quadLine(jj).x(2).XData = squeeze(obj.quadData(jj).propeller(:,1,2));
                obj.quadLine(jj).x(2).YData = squeeze(obj.quadData(jj).propeller(:,2,2));
                obj.quadLine(jj).x(3).XData = squeeze(obj.quadData(jj).propeller(:,1,3));
                obj.quadLine(jj).x(3).YData = squeeze(obj.quadData(jj).propeller(:,2,3));
                obj.quadLine(jj).x(4).XData = squeeze(obj.quadData(jj).propeller(:,1,4));
                obj.quadLine(jj).x(4).YData = squeeze(obj.quadData(jj).propeller(:,2,4));

                obj.quadLine(jj).x(5).XData = squeeze(obj.quadData(jj).frame(:,1,1));
                obj.quadLine(jj).x(5).YData = squeeze(obj.quadData(jj).frame(:,2,1));
                obj.quadLine(jj).x(6).XData = squeeze(obj.quadData(jj).frame(:,1,2));
                obj.quadLine(jj).x(6).YData = squeeze(obj.quadData(jj).frame(:,2,2));
                obj.quadLine(jj).path.XData = [obj.quadLine(jj).path.XData xpath(ii,jj)];
                obj.quadLine(jj).path.YData = [obj.quadLine(jj).path.YData ypath(ii,jj)];

                obj.quadLine(jj).orientation.XData = obj.quadData(jj).orient(:,1);
                obj.quadLine(jj).orientation.YData = obj.quadData(jj).orient(:,2);
                yawangle.XData = [yawangle.XData dt*ii];
                yawangle.YData = [yawangle.YData psipath(ii)];
                massLine.XData = [massLine.XData dt*ii];
                massLine.YData = [massLine.YData mass(ii)];
                end
                % Update the wind field speed
                U = 1+0.5*rand(size(X));
                V = -1+0.5*rand(size(Y));
                g.UData = U;
                g.VData = V;
                t2 = toc(t1);
                if t2<0.04
                    pause(0.04-t2)
                else 
                    pause(0.01);
                end
                t1 = tic;
            end
            disp('Done')
        end
    end
end

