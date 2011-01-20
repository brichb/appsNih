#ifndef MY_QFILE_DIALOG_H
#define MY_QFILE_DIALOG_H

#include <qfiledialog.h>
#include <qpushbutton.h>
#include <qlabel.h>

class myQFileDialog : public QFileDialog 
{
	Q_OBJECT

	public:
		QStringList selectedDirs;
		QLabel* label;
		myQFileDialog( QWidget * parent = 0, const char * name = 0, bool modal = FALSE )
			: QFileDialog(parent, name, modal)
		  {
		       label = new QLabel( "dirs", this );
		       QPushButton* b = new QPushButton( this, "GetDir" );
				b->setText("Add Dir");
				addWidgets( label, NULL, b );
				connect(b, SIGNAL(clicked()), this, SLOT(getDirs()) );
				this->setMode(QFileDialog::DirectoryOnly );
		};
public slots:
		  void getDirs()
		  {
			  QString sf = selectedFile();
			  selectedDirs.push_back(sf);
			  QString allDirs = selectedDirs.join("\n");
			  label->setText(allDirs);

		  };
};


#endif
