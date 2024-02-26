import { BehaviorSubject, firstValueFrom, Observable } from 'rxjs';
import { DomSanitizer } from '@angular/platform-browser';
import { HttpClient } from '@angular/common/http';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { Output, RawOutput } from './output';
import { OutputServiceInterface } from './output.service.interface';
import { UserIdService } from './user-id.service';

@Injectable({
  providedIn: 'root',
})
export class OutputService implements OutputServiceInterface {
  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  private userId: string | null = null;

  constructor(private http: HttpClient, private userIdService: UserIdService, private snackBar: MatSnackBar, private sanitizer: DomSanitizer) {
    this.userIdService.userId.subscribe(
      (res) => { this.userId = res; },
    );
  }

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  async executeBlock(block: Block): Promise<void> {
    if (this.userId === null) {
      this.snackBar.open('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    } else {
      const url = 'http://127.0.0.1:5000/' + block.blockId;
      type Params = { [key: string]: string | number };
      const params: Params = {};
      params['user_id'] = this.userId;
      for (let i = 0; i < block.parameters.length; i++) {
        params[block.parameters[i].key] = block.parameters[i].value;
      }
      try {
        const outputResponse: RawOutput = await firstValueFrom(this.http.get<RawOutput>(url, { params: params }));
        let processedResponse: Output = {};
        if (outputResponse.text) {
          processedResponse = outputResponse;
        } else if (outputResponse.img) {
          const imageString = outputResponse.img as string;
          const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
          const objectURL = 'data:image/png;base64,' + processedString;
          const newImg = this.sanitizer.bypassSecurityTrustUrl(objectURL);
          processedResponse = { 
            img: newImg,
            alttext: outputResponse.alttext,
          };
        }
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } catch (e) {
        this.snackBar.open('Not a valid order, please check the blocks and try again', 'Close', { duration: 5000 });
      }
    }
  }
}
