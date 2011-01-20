#include "NIH_ActiveLocationModel.h"

#include <qstringlist.h>
#include "XmlReader.h"

ObjectLocationModel::ObjectLocationModel()
{
	center = Vec3(0,0,0);
	size = Vec3(0,0,0);
	orientation = UnitRotation;
	state = 0;
	shape = 0;

	strcpy(name, "NULL");
}

ObjectLocationModel::~ObjectLocationModel()
{
}

RelativeLocationModel::RelativeLocationModel()
{
	relativeLocation = Vec3(0,0,0);
	relativeOrientation = UnitRotation;
	state = 0;
	strcpy(name, "NULL");

	src = target = NULL;
}

RelativeLocationModel::~RelativeLocationModel()
{
}

VertebraLocationModel::VertebraLocationModel()
: ObjectLocationModel()
{
}

VertebraLocationModel::~VertebraLocationModel()
{
}

OrganLocationModel::OrganLocationModel()
: ObjectLocationModel()
{
///	toL1 = toL2 = toT11 = toT12 = NULL;
	for(int i=0; i<6; i++)
	{
		toTL[i].relativeLocation = Vec3(0,0,0);
	}
}

OrganLocationModel::~OrganLocationModel()
{
///	if(toL1!=NULL) delete toL1;
///	if(toL2!=NULL) delete toL1;
///	if(toT12!=NULL) delete toT12;
///	if(toT11!=NULL) delete toT11;
}

LocationFramework::LocationFramework()
{
	int i;
	for(i=0; i<17; i++) tL[i]=NULL;
	for(i=0; i<5; i++) organ[i]=NULL;

}

LocationFramework::~LocationFramework()
{	
	int i;
	for(i=0; i<17; i++) 
		if(tL[i]!=NULL) delete tL[i];

	for(i=0; i<5; i++) 
		if(organ[i]!=NULL) delete organ[i];

}

ActiveLocationModel_PDM::ActiveLocationModel_PDM()
{
	numTrainingVec = 0;
	vecLength = 0;
	pdmDegree = 0;
	strcpy(referenceModelFileName, "");
}

ActiveLocationModel_PDM::~ActiveLocationModel_PDM()
{
}

int ActiveLocationModel_PDM::AddTrainingVec(CIS_Vector_JY_double &tVec)
{
	if(numTrainingVec==0) vecLength = tVec.GetSize();
	else if(vecLength!=tVec.GetSize())
	{
		// all training vectors should have same length
		return 0;
	}
	trainingVec[numTrainingVec] = tVec;
	numTrainingVec ++;
	return numTrainingVec;
}

int ActiveLocationModel_PDM::TrainPDM()
{
	// should have at least two training Vectors
	if(numTrainingVec<2) return 1;

	CIS_Vector_JY_double Y_d;
	CIS_Matrix_JY_double S, D, Dt;
	int i, j;

	meanVec.SetSize(vecLength);
	orgVec.SetSize(vecLength);

	// compute average vector
	meanVec = trainingVec[0];
	for(i=1; i<numTrainingVec; i++) meanVec += trainingVec[i];
	meanVec /= (double)numTrainingVec;

	// set Y_org
	orgVec = trainingVec[0];

	S.SetSize(numTrainingVec, numTrainingVec);
	D.SetSize(vecLength, numTrainingVec);
	Dt.SetSize(numTrainingVec, vecLength);
	Y_d.SetSize(vecLength);

	// build the covariance matrix
	for(i=0; i<numTrainingVec; i++)
	{
		Y_d = trainingVec[i];
		Y_d -= meanVec;

		for(j=0; j<vecLength; j++)
		{
			D(j, i) = Y_d[j];
		}
	}

	Dt = D.Transpose();

	S = Dt * D;

	S /= (double)(numTrainingVec);

	// compute the eigen vector
	int nrot;
	CIS_Matrix_JY_double e_vector;
	CIS_Vector_JY_double e_value;

	e_value.SetSize(numTrainingVec);
	e_vector.SetSize(numTrainingVec, numTrainingVec);
	NR_Eigen_jacobi(S, numTrainingVec, e_value, e_vector, &nrot);

	// compute the real e-vector;
	pdmModel.SetSize(vecLength, numTrainingVec);
	pdmModel = D * e_vector;

	// normalize the P vector
	double norm;
	for(i=0; i<numTrainingVec; i++)
	{
		norm=0;
		for(j=0; j<vecLength; j++)
		{
			norm += pdmModel(j,i)*pdmModel(j,i);
		}

		if(norm!=0)
		{
			norm = sqrt(norm);
			for(j=0; j<vecLength; j++)
			{
				pdmModel(j, i) /= norm;
			}	 // for j
		}	// if norm
	}	// for i

	pdmDegree = numTrainingVec;
	pdmRange.SetSize(numTrainingVec);
	pdmValue.SetSize(numTrainingVec);
	// 3 standard deviation of PCA model
	for(i=0 ;i<numTrainingVec; i++) 
	{
		if(e_value[i]<1e-3) pdmRange[i]=0;
		else pdmRange[i] = sqrt(e_value[i]);
		pdmValue[i] = 0;
	}

	return 0;
}

// intanstiate an PDMmodel
CIS_Vector_JY_double ActiveLocationModel_PDM::InstantiatePDM(doubleDynArray &b)
{
	CIS_Vector_JY_double ret;

	if(b.GetSize()>pdmModel.Num_Cols())		
	{
		printf("Error: too many parameters for instantiation\n");
		return ret;
	}

	CIS_Matrix_JY_double P_;
	int i, j;
	CIS_Vector_JY_double b_;

	int degree=b.GetSize();
	b_.SetSize(degree);
	P_.SetSize(pdmModel.Num_Rows(), degree);

	for(i=0; i<degree; i++) 
	{
		b_[i] = b[i];
		for(j=0; j<P_.Num_Rows(); j++)
		{
			P_(j, i) = pdmModel(j,i);
		}
	}

	ret = meanVec + P_*b_;

	return ret;
}

CIS_Vector_JY_double ActiveLocationModel_PDM::InstantiatePDM()
{
	return InstantiatePDM(pdmValue);
}

// given the feature vector, computer the pdmPara meters
int ActiveLocationModel_PDM::InversePDM(CIS_Vector_JY_double &tVec, doubleDynArray &b)
{
	CIS_Matrix_JY_double P_, P_t, P_square;
	int i, j;
	CIS_Vector_JY_double b_;

	int degree=pdmDegree;
	b_.SetSize(degree);
	P_.SetSize(pdmModel.Num_Rows(), degree);

	for(i=0; i<degree; i++) 
	{
		for(j=0; j<P_.Num_Rows(); j++)
		{
			P_(j, i) = pdmModel(j,i);
		}
	}

	P_t = P_.Transpose();

	P_square = P_t * P_;
	P_square = P_square.Inverse();

	b_ = P_square * P_t * (tVec-meanVec);

	b.SetSize(degree);

	for(i=0; i<degree; i++) b[i]=b_[i];
	return 0;
}

int ActiveLocationModel_PDM::SaveALMModel(const char *fn)
{
	int i, j;
	FILE *fp = fopen(fn,"w");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open xml model file %s\n",fn);
		return 1;
	}

	fprintf(fp, "<ALMModel>\n");

	// Save information about the reference model
	fprintf(fp,"\t<referenceModel filename=\"%s\">\n", referenceModelFileName);
	fprintf(fp,"\t\t<referenceModelSize> %f %f %f </referenceModelSize>\n", referenceModelSize.x, referenceModelSize.y, referenceModelSize.z);
	// save vertex and organ list
	fprintf(fp,"\t\t<vertebraIndexList size=\"%d\">\n", vertebraIndexList.GetSize());
	fprintf(fp,"\t\t\t<v> ");
	for(i=0; i<vertebraIndexList.GetSize(); i++)
	{
		fprintf(fp,"%d ", vertebraIndexList[i]);
	}
	fprintf(fp,"</v>\n");
	fprintf(fp,"\t\t</vertebraIndexList>\n");
	fprintf(fp,"\t\t<organIndexList size=\"%d\">\n", organIndexList.GetSize());
	fprintf(fp,"\t\t\t<v> ");
	for(i=0; i<organIndexList.GetSize(); i++)
	{
		fprintf(fp,"%d ", organIndexList[i]);
	}
	fprintf(fp,"</v>\n");
	fprintf(fp,"\t\t</organIndexList>\n");
	fprintf(fp,"\t</referenceModel>\n");

	// Save the pdmModel
	fprintf(fp,"\t<pdmModel vecLength=\"%d\" pdmDegree=\"%d\">\n", vecLength, pdmDegree);

	fprintf(fp,"\t\t<meanVec>\n");
	for(i=0; i<meanVec.GetSize(); i++)
	{
		fprintf(fp,"\t\t\t<v> %f </v>\n", meanVec[i]);
	}
	fprintf(fp,"\t\t</meanVec>\n");

	fprintf(fp,"\t\t<orgVec>\n");
	for(i=0; i<orgVec.GetSize(); i++)
	{
		fprintf(fp,"\t\t\t<v> %f </v>\n", orgVec[i]);
	}
	fprintf(fp,"\t\t</orgVec>\n");

	fprintf(fp,"\t\t<pdmRange>\n");
	for(i=0; i<pdmRange.GetSize(); i++)
	{
		fprintf(fp,"\t\t\t<v> %f </v>\n", pdmRange[i]);
	}
	fprintf(fp,"\t\t</pdmRange>\n");

	fprintf(fp,"\t\t<pdmValue>\n");
	for(i=0; i<pdmModel.Num_Cols(); i++)
	{
		fprintf(fp,"\t\t\t<v>");
		for(j=0; j<pdmModel.Num_Rows(); j++) fprintf(fp, "%f ", pdmModel(j, i));
		fprintf(fp,"</v>\n");
	}
	fprintf(fp,"\t\t</pdmValue>\n");

	fprintf(fp,"\t</pdmModel>\n");

	fprintf(fp, "</ALMModel>\n");
	fclose(fp);

	return 0;
}

int ActiveLocationModel_PDM::LoadALMModel(const char *fn)
{
	FILE *fp = fopen(fn,"r");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open xml model file %s\n",fn);
		return -1;
	}
	fclose(fp);

	int i, j;
	QStringList valueList;
    // read the xml file
    XmlReader xmlReader(fn);
    QDomElement  ALMModelElement = xmlReader.getDomElementUnique("ALMModel", xmlReader.getTopNode());

	QDomElement referenceModelElement = xmlReader.getDomElementUnique("referenceModel",ALMModelElement);
	strcpy(referenceModelFileName, referenceModelElement.attribute("filename").ascii());
	valueList = QStringList::split(' ', referenceModelElement.elementsByTagName("referenceModelSize").item(0).toElement().text().simplifyWhiteSpace());
	referenceModelSize.x = valueList[0].toDouble();
	referenceModelSize.y = valueList[1].toDouble();
	referenceModelSize.z = valueList[2].toDouble();

	// get the vertebraIndexList
	QDomElement vertebraIndexListElement = xmlReader.getDomElementUnique("vertebraIndexList",referenceModelElement);
	vertebraIndexList.SetSize(vertebraIndexListElement.attribute("size").toInt());
	valueList = QStringList::split(' ', vertebraIndexListElement.elementsByTagName("v").item(0).toElement().text().simplifyWhiteSpace());
	for(i=0; i<valueList.size(); i++)
	{
		vertebraIndexList[i] = valueList[i].toInt();
	}

	// get the organIndexList
	QDomElement organIndexListElement = xmlReader.getDomElementUnique("organIndexList",referenceModelElement);
	organIndexList.SetSize(organIndexListElement.attribute("size").toInt());
	valueList = QStringList::split(' ', organIndexListElement.elementsByTagName("v").item(0).toElement().text().simplifyWhiteSpace());
	for(i=0; i<valueList.size(); i++)
	{
		organIndexList[i] = valueList[i].toInt();
	}

	QDomElement pdmModelElement = xmlReader.getDomElementUnique("pdmModel",ALMModelElement);
	vecLength = pdmModelElement.attribute("vecLength").toInt();
	pdmDegree = pdmModelElement.attribute("pdmDegree").toInt();

	meanVec.SetSize(vecLength);
	QDomElement meanVecElement = xmlReader.getDomElementUnique("meanVec",pdmModelElement);
    QDomNodeList vElementList    = xmlReader.getDomNodeList("v", meanVecElement);
	for(i=0; i<vElementList.count(); i++)
	{
		meanVec[i] = vElementList.item(i).toElement().text().toDouble();
	}

	orgVec.SetSize(vecLength);
	QDomElement orgVecElement = xmlReader.getDomElementUnique("orgVec",pdmModelElement);
    vElementList    = xmlReader.getDomNodeList("v", orgVecElement);
	for(i=0; i<vElementList.count(); i++)
	{
		orgVec[i] = vElementList.item(i).toElement().text().toDouble();
	}
	
	pdmRange.SetSize(pdmDegree);
	pdmValue.SetSize(pdmDegree);
	QDomElement pdmRangeElement = xmlReader.getDomElementUnique("pdmRange",pdmModelElement);
    vElementList    = xmlReader.getDomNodeList("v", pdmRangeElement);
	for(i=0; i<vElementList.count(); i++)
	{
		pdmRange[i] = vElementList.item(i).toElement().text().toDouble();
	}

	pdmModel.SetSize(vecLength, pdmDegree);
	QDomElement pdmValueElement = xmlReader.getDomElementUnique("pdmValue",pdmModelElement);
    vElementList    = xmlReader.getDomNodeList("v", pdmValueElement);
	for(i=0; i<vElementList.count(); i++)
	{
		valueList = QStringList::split(' ', vElementList.item(i).toElement().text().simplifyWhiteSpace());
		for(j=0; j<valueList.size(); j++) pdmModel(j, i)=valueList[j].toDouble();
	}

	
	return 0;
}