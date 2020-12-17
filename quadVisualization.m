classdef quadVisualization
    %quadVisualization Visualize the quadcopter based on its dimentsion
    %   Detailed explanation goes here
    
    properties
        L = 0.4*sqrt(2)
        r = 0.2
        d = 0.4*sqrt(2)/10
        xO = 0 % center of the quadcopter
        yO = 0
        psi = 0 % course of the quadcopter
        propeller
        frame
        orient
    end
    
    methods
        function obj = quadVisualization(L,r,d)
            %quadVisualization Construct an instance of this class
            %   Detailed explanation goes here
            switch(nargin)
                case 1
                    obj.L = L;
                case 2
                    obj.L = L;
                    obj.r = r;
                case 3
                    obj.L = L;
                    obj.r = r;
                    obj.d = d;
            end
            xOprop = obj.L*cos(pi/4+obj.psi);
            yOprop = obj.L*sin(pi/4+obj.psi);
            
            phi= 0:0.1:2*pi+0.1;
            xCprop = obj.r*cos(phi);
            yCprop = obj.r*sin(phi);
            
            obj.propeller(:,:,1) = [obj.xO+xOprop+xCprop;obj.yO+yOprop+yCprop]';
            obj.propeller(:,:,2) = [obj.xO-yOprop+xCprop;obj.yO+xOprop+yCprop]';
            obj.propeller(:,:,3) = [obj.xO-xOprop+xCprop;obj.yO-yOprop+yCprop]';
            obj.propeller(:,:,4) = [obj.xO+yOprop+xCprop;obj.yO-xOprop+yCprop]';
           
            xFrame = [obj.xO+(-xOprop:2*xOprop/100:xOprop)+obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-xOprop:2*xOprop/100:xOprop)-obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-yOprop:2*yOprop/100:yOprop)-obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-yOprop:2*yOprop/100:yOprop)+obj.d*sin(pi/4+obj.psi)];
            yFrame = [obj.yO+(-yOprop:2*yOprop/100:yOprop)-obj.d*cos(pi/4+obj.psi);...
                obj.yO+(-yOprop:2*yOprop/100:yOprop)+obj.d*cos(pi/4+obj.psi);...
                obj.yO+(xOprop:-2*xOprop/100:-xOprop)-obj.d*cos(pi/4+obj.psi);...
                obj.yO+(xOprop:-2*xOprop/100:-xOprop)+obj.d*cos(pi/4+obj.psi)];
            obj.frame(:,:,1) = [[xFrame(1,:) xFrame(2,end:-1:1) xFrame(1,1)]; [yFrame(1,:) yFrame(2,end:-1:1) yFrame(1,1)]]';
            obj.frame(:,:,2) = [[xFrame(3,:) xFrame(4,end:-1:1) xFrame(3,1)]; [yFrame(3,:) yFrame(4,end:-1:1) yFrame(3,1)]]';
            obj.orient = [[obj.xO obj.xO+20*obj.d*cos(obj.psi)]; [obj.yO obj.yO+20*obj.d*sin(obj.psi)]]';
        end
        
        function obj = update(obj, xO, yO, psi)
            %update update propeller, frame and orientation
            obj.xO = xO;
            obj.yO = yO;
            obj.psi = psi;
            xOprop = obj.L*cos(pi/4+obj.psi);
            yOprop = obj.L*sin(pi/4+obj.psi);
            
            phi= 0:0.1:2*pi+0.1;
            xCprop = obj.r*cos(phi);
            yCprop = obj.r*sin(phi);
            
            obj.propeller(:,:,1) = [obj.xO+xOprop+xCprop;obj.yO+yOprop+yCprop]';
            obj.propeller(:,:,2) = [obj.xO-yOprop+xCprop;obj.yO+xOprop+yCprop]';
            obj.propeller(:,:,3) = [obj.xO-xOprop+xCprop;obj.yO-yOprop+yCprop]';
            obj.propeller(:,:,4) = [obj.xO+yOprop+xCprop;obj.yO-xOprop+yCprop]';
           
            xFrame = [obj.xO+(-xOprop:2*xOprop/100:xOprop)+obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-xOprop:2*xOprop/100:xOprop)-obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-yOprop:2*yOprop/100:yOprop)-obj.d*sin(pi/4+obj.psi);...
                obj.xO+(-yOprop:2*yOprop/100:yOprop)+obj.d*sin(pi/4+obj.psi)];
            yFrame = [obj.yO+(-yOprop:2*yOprop/100:yOprop)-obj.d*cos(pi/4+obj.psi);...
                obj.yO+(-yOprop:2*yOprop/100:yOprop)+obj.d*cos(pi/4+obj.psi);...
                obj.yO+(xOprop:-2*xOprop/100:-xOprop)-obj.d*cos(pi/4+obj.psi);...
                obj.yO+(xOprop:-2*xOprop/100:-xOprop)+obj.d*cos(pi/4+obj.psi)];
            obj.frame(:,:,1) = [[xFrame(1,:) xFrame(2,end:-1:1) xFrame(1,1)]; [yFrame(1,:) yFrame(2,end:-1:1) yFrame(1,1)]]';
            obj.frame(:,:,2) = [[xFrame(3,:) xFrame(4,end:-1:1) xFrame(3,1)]; [yFrame(3,:) yFrame(4,end:-1:1) yFrame(3,1)]]';
            obj.orient = [[obj.xO obj.xO+20*obj.d*cos(obj.psi)]; [obj.yO obj.yO+20*obj.d*sin(obj.psi)]]';
        end
        
        
    end
end

