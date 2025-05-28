function iou = boxIoU(l1, t1, w1, h1, l2, t2, w2, h2)
    area1 = w1 * h1;
	area2 = w2 * h2;

	x_overlap = max(0, min(l1 + w1, l2 + w2) - max(l1, l2));
	y_overlap = max(0, min(t1 + h1, t2 + h2) - max(t1, t2));
	intersectionArea = x_overlap*y_overlap;
	unionArea = area1 + area2 - intersectionArea;
    iou = intersectionArea / unionArea;
end