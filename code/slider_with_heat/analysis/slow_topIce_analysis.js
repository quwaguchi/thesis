const fs = require('fs');
const path = require('path');

async function processFile(fileNum, profFolderPath, outputFolderPath) {
   //ファイル名は output_0000.prof ~ output_fileNum.prof
   //profFolderPathは上記profファイルを格納しているフォルダーの絶対パス

   const eps = 1e-10;
   let outputStr = `Prof Folder Path: ${profFolderPath}\n`;

   for (let i = 0; i <= fileNum; i++) {
      const fileIndex = ('0000' + i).slice(-4);
      const filePath = `${profFolderPath}/output_${fileIndex}.prof`;

      const data = await fs.readFileSync(filePath, 'utf8');
      const lines = data.split('\n');
      const time = lines[0];

      let tempArr = [];

      for (let j=2; j<lines.length; j++){
         const particleData = lines[j].split(' ');
         const x = parseFloat(particleData[1]);
         const y = parseFloat(particleData[2]);
         const temp = parseFloat(particleData[5]);

         if (x<0.000000-eps || 1.500000+eps<=x) continue;

         if (Math.abs(y-0.014000)<eps){
            tempArr.push(temp);
         }
      }

      const tempStr = tempArr.join(', ');
      outputStr = outputStr.concat(`${time}, ${tempArr}\n`);
   }

   const outputFilePath = `${outputFolderPath}/topIce_analysis.csv`;
   fs.writeFileSync(outputFilePath, outputStr);
}


// processFile(fileNum, profFolderPath, outputFolderPath)
processFile(
   1049,
   "/Users/quwaguchi/thesis/thesis/code/slider_with_heat/with_friction/(-1,-1)/slow/prof",
   "/Users/quwaguchi/thesis/thesis/code/slider_with_heat/with_friction/(-1,-1)/slow"
)
