import { Component } from '@angular/core';
import { DomSanitizer } from '@angular/platform-browser';

import { Output } from './../output';
import { OutputService } from './../output.service';

@Component({
  selector: 'app-output-display',
  templateUrl: './output-display.component.html',
  styleUrls: ['./output-display.component.css'],
})

export class OutputDisplayComponent {
  outputList: Output[] = [];

  constructor(private outputService: OutputService, sanitizer: DomSanitizer) {
    this.outputService.outputs.subscribe(
      (res) => { 
        this.outputList = res.map(out => {
          if (out.text) {
            return out;
          } else if (out.img) {
            const imageString = out.img as string;
            let fixedData = imageString.substring(2, imageString.length-3);
            fixedData = fixedData.replace(/\\n/g, '');
            const objectURL = 'data:image/png;base64,' + fixedData;
            const newImg = sanitizer.bypassSecurityTrustUrl(objectURL);
            return {
              img: newImg
            }; 
          } else {
            return out;
          }
        }); 
      }
    );
  }
}
