const fs = require('fs');
const path = require('path');

async function processFile(fileNum, profFolderPath, outputFolderPath) {
   //ファイル名は output_0000.prof ~ output_fileNum.prof
   //profFolderPathは上記profファイルを格納しているフォルダーの絶対パス

   const eps = 1e-10;
   let output = `Prof Folder Path: ${profFolderPath}\n`;

   for (let i = 0; i <= fileNum; i++) {
      const fileIndex = ('0000' + i).slice(-4);
      const filePath = `${profFolderPath}/output_${fileIndex}.prof`;

      const data = await fs.readFileSync(filePath, 'utf8');
      const lines = data.split('\n');
      const time = lines[0];

      let L1 = 0;
      let L1_FnS = 0;
      let L2 = 0;
      let L2_FnS = 0;
      let L3 = 0;
      let L3_FnS = 0;

      for (let j=2; j<lines.length; j++){
         const particleData = lines[j].split(' ');
         const particleType = parseInt(particleData[0]);
         const x = parseFloat(particleData[1]);
         const y = parseFloat(particleData[2]);

         if (x<0.000000-eps || 1.500000+eps<=x) continue;

         if (Math.abs(y-0.014000)<eps){
            L1++;
            if (particleType===1){
               L1_FnS++;
            }
         }else if (Math.abs(y-0.013000)<eps){
            L2++;
            if (particleType===1){
               L2_FnS++;
            }
         }else if (Math.abs(y-0.012000)<eps){
            L3++;
            if (particleType===1){
               L3_FnS++;
            }
         }
      }

      output = output.concat(`${time}, ${L1_FnS/L1}, ${L2_FnS/L2}, ${L3_FnS/L3}\n`);
   }

   const outputFilePath = `${outputFolderPath}/FnS_analysis.csv`;
   fs.writeFileSync(outputFilePath, output);
}


// processFile(fileNum, profFolderPath, outputFolderPath)
processFile(
   1049,
   "/Users/quwaguchi/thesis/thesis/code/slider_with_heat/with_friction/(-1,-1)/slow/prof",
   "/Users/quwaguchi/thesis/thesis/code/slider_with_heat/with_friction/(-1,-1)/slow"
)
