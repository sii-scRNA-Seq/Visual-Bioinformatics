import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { BehaviorSubject, Observable, firstValueFrom } from 'rxjs';
import { Output } from './output';
import { Block } from './block.service';
import { MatSnackBar } from '@angular/material/snack-bar';

@Injectable({
  providedIn: 'root',
})

export class OutputService {
  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  constructor(private http: HttpClient, private snackBar: MatSnackBar) { }

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  async executeBlock(block: Block): Promise<void> {
    const url = 'http://127.0.0.1:5000/' + block.blockId;
    type Params = {[key: string] : number};
    const params: Params = {};
    for (let i=0; i < block.parameters.length; i++) {
      params[ block.parameters[i].key ] = block.parameters[i].value;
    }
    let outputResponse: Output = {text: '', other: ''};
    try {
      outputResponse = await firstValueFrom(this.http.get<Output>(url, {params: params}));
    } catch(e) {
      outputResponse = {text: '', other: ''};
      this.snackBar.open('Not a valid order, please check the blocks and try again', 'Close', { duration: 5000 });
    } finally {
      const outputs = this.outputs$.getValue();
      outputs.push(outputResponse);
      this.outputs$.next(outputs);
    }
  }
}
