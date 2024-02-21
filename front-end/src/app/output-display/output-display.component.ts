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
        this.outputList = res.map(output => {
          if (output.text) {
            return output;
          } else if (output.img) {
            const imageString = output.img as string;
            const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
            const objectURL = 'data:image/png;base64,' + processedString;
            const newImg = sanitizer.bypassSecurityTrustUrl(objectURL);
            return {
              img: newImg
            }; 
          } else {
            return output;
          }
        }); 
      }
    );
  }
}
