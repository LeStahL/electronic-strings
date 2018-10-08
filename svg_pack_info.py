pack svg into texture
layout: 
- 2 Byte, short: Number of contained elements N
- N Blocks consisting of element data, consisting of
	- 1 Byte, char: Element ID.
		-0: Line
		-1: Spline
		-2: Rectangle
		-3: Circle
		-4: Polygon
		-5: Stroke
		-6: Text
	- Element control data.
		-0: Line
			-8 Bytes, float16[4], {x1,y1,x2,y2}
			-4 Bytes, float8[4], color {r,g,b,a}
			-2 Bytes, float16, width
		-1: Spline
			-12 Bytes, float16[6], {x1,y1,x2,y2,x3,y3}
			-4 Bytes, float8[4], color
			-2 Bytes, float16, width
		-2: Rectangle
			-4 Bytes, float16[2], {x1,y1}
			-4 Bytes, float16[2], {width, height}

